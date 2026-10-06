"""Aşama 3 — koşu sonrası analiz: A2/A0, enerji, Sigma(R), disk kalınlığı, animasyon.

Yalnızca runs/<ad>/ altındaki dosyaları okur; simülasyon koşmaz, koşu sırasında
çağrılmaz. Zaman serileri (E, A2/A0, kütle merkezi) koşu sırasında
timeseries.npz'ye yazılmıştır; snapshot'lar yalnızca uzaysal ölçümler için okunur.

Birimler: konum : kpc, hız : kpc/Myr, zaman : Myr, kütle : M_sun,
Sigma : M_sun/kpc^2, enerji : M_sun kpc^2/Myr^2.

Not: config.py bu aşamada değiştirilmediği için analiz ayarları aşağıda
modül sabitleri olarak durur (CLAUDE.md §8'e göre ileride config.py'ye taşınmalı).
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .config import A2_R_MAX
from .nbody import bar_strength
from .units import KMS_KPMY

# Analiz ayarları (kpc, Myr, boyutsuz)
SCALE_FIT_RANGE = (3.0, 10.0)   # kpc, üstel Sigma uyumunun yarıçap aralığı
SCALE_FIT_BINS = 14             # uyum aralığındaki halka sayısı
PROFILE_R_MAX = 20.0            # kpc, Sigma(R) profilinin dış sınırı
PROFILE_BINS = 40
BAR_THRESHOLD = 0.2             # A2/A0 bu değeri aşınca "çubuk var" denir
BAR_SMOOTH = 100.0              # Myr, çubuk zamanı için hareketli ortalama penceresi
CLUMP_PIXEL = 0.25              # kpc, topaklanma ölçümünde piksel boyu
CLUMP_CONTRAST = 3.0            # azimut ortalamasının bu katından yoğun piksel = topak
CLUMP_R_MIN = 2.0               # kpc, merkez (çubuk/şişkinlik) dışarıda tutulur
MAP_HALF_WIDTH = 15.0           # kpc, yüzey yoğunluğu haritalarının yarı genişliği


# ==========================================================================
# Koşu okuma
# ==========================================================================

@dataclass
class Run:
    """Bir koşu klasörü. t_snap : (S,) Myr; ts : zaman serileri; config : dict."""
    name: str
    path: Path
    t_snap: np.ndarray
    mass: np.ndarray
    ts: dict
    config: dict

    @property
    def n_snap(self):
        return len(self.t_snap)

    def snapshot(self, k):
        """k. snapshot: (N,6) float64 [x,y,z,vx,vy,vz] — kpc, kpc/Myr."""
        return np.load(self.path / f"snap_{k:04d}.npy").astype(np.float64)

    def index_at(self, t):
        """t [Myr]'ye en yakın snapshot indeksi."""
        return int(np.argmin(np.abs(self.t_snap - t)))


def is_run_dir(path):
    """Aşama 3 biçiminde (timeseries.npz + config.json) bir koşu klasörü mü?"""
    p = Path(path)
    return (p / "timeseries.npz").exists() and (p / "config.json").exists() and (p / "times.npy").exists()


def is_complete(path):
    """Koşu bitmiş mi: snapshot sayısı t_end / snap_interval + 1'e ulaşmış olmalı."""
    p = Path(path)
    if not is_run_dir(p):
        return False
    cfg = json.loads((p / "config.json").read_text())
    expected = int(round(cfg["t_end"] / cfg["numerical"]["snap_interval"])) + 1
    return len(list(p.glob("snap_*.npy"))) >= expected


def list_runs(runs_dir="runs"):
    """runs_dir altındaki tamamlanmış Aşama 3 biçimli koşuların adları (sıralı)."""
    root = Path(runs_dir)
    if not root.exists():
        return []
    return sorted(p.name for p in root.iterdir() if p.is_dir() and is_complete(p))


def load_run(name, runs_dir="runs"):
    p = Path(runs_dir) / name
    return Run(
        name=name,
        path=p,
        t_snap=np.load(p / "times.npy"),
        mass=np.load(p / "mass.npy"),
        ts=dict(np.load(p / "timeseries.npz")),
        config=json.loads((p / "config.json").read_text()),
    )


# ==========================================================================
# Zaman serileri
# ==========================================================================

def energy_error(run: Run):
    """(t [Myr], (E - E0)/|E0|)."""
    E = run.ts["E"]
    return run.ts["t"], (E - E[0]) / abs(E[0])


def energy_drift_rate(run: Run, t_min=200.0):
    """t >= t_min aralığında (E-E0)/|E0| doğrusal eğimi, 1/Myr."""
    t, dE = energy_error(run)
    s = t >= t_min
    return float(np.polyfit(t[s], dE[s], 1)[0])


def bar_series(run: Run):
    """(t [Myr], A2/A0) — koşu sırasında R < A2_R_MAX içinde ölçülen."""
    return run.ts["t"], run.ts["A2"]


def running_mean(t, y, window):
    """Eşit aralıklı t üzerinde, window [Myr] genişliğinde merkezli hareketli ortalama."""
    dt = np.median(np.diff(t))
    n = max(1, int(round(window / dt)))
    kernel = np.ones(n) / n
    pad = np.pad(y, (n // 2, n - 1 - n // 2), mode="edge")
    return np.convolve(pad, kernel, mode="valid")


def bar_formation_time(t, A2, threshold=BAR_THRESHOLD, window=BAR_SMOOTH):
    """Hareketli ortalaması (window Myr) threshold'u ilk aştığı an [Myr]; yoksa None."""
    sm = running_mean(t, A2, window)
    above = np.flatnonzero(sm > threshold)
    return float(t[above[0]]) if above.size else None


def com_drift(run: Run):
    """Kütle merkezi kayması: (max |r_com| [kpc], max |v_com| [kpc/Myr])."""
    return (float(np.linalg.norm(run.ts["com"], axis=1).max()),
            float(np.linalg.norm(run.ts["vcom"], axis=1).max()))


# ==========================================================================
# Snapshot ölçümleri
# ==========================================================================

def surface_density_profile(snap, mass, R_max=PROFILE_R_MAX, n_bins=PROFILE_BINS):
    """Halka bazında Sigma(R). dönüş : (R_mid [kpc], Sigma [M_sun/kpc^2])."""
    R = np.hypot(snap[:, 0], snap[:, 1])
    e = np.linspace(0.0, R_max, n_bins + 1)
    S = np.histogram(R, e, weights=mass)[0] / (np.pi * (e[1:] ** 2 - e[:-1] ** 2))
    return 0.5 * (e[1:] + e[:-1]), S


def scale_length(snap, mass, R_range=SCALE_FIT_RANGE, n_bins=SCALE_FIT_BINS):
    """Sigma(R)'ye R_range içinde üstel uyumdan R_d [kpc]."""
    R = np.hypot(snap[:, 0], snap[:, 1])
    e = np.linspace(*R_range, n_bins + 1)
    S = np.histogram(R, e, weights=mass)[0] / (np.pi * (e[1:] ** 2 - e[:-1] ** 2))
    ok = S > 0
    return float(-1.0 / np.polyfit(0.5 * (e[1:] + e[:-1])[ok], np.log(S[ok]), 1)[0])


def disk_thickness(snap):
    """Disk kalınlığı: medyan |z - <z>| [kpc] (sech^2 için 0.549 z_0)."""
    z = snap[:, 2]
    return float(np.median(np.abs(z - z.mean())))


def bar_strength_snapshot(snap, mass, R_max=A2_R_MAX):
    """Snapshot'tan A2/A0 (nbody.bar_strength ile aynı tanım)."""
    return bar_strength(snap[:, :3], mass, R_max)


def bar_strength_profile(snap, mass, R_max=PROFILE_R_MAX / 2, n_bins=PROFILE_BINS // 2):
    """Halka bazında A2/A0(R). dönüş : (R_mid [kpc], A2/A0)."""
    x, y = snap[:, 0], snap[:, 1]
    R, ph = np.hypot(x, y), np.arctan2(y, x)
    e = np.linspace(0.0, R_max, n_bins + 1)
    idx = np.digitize(R, e) - 1
    ok = (idx >= 0) & (idx < n_bins)
    a0 = np.bincount(idx[ok], weights=mass[ok], minlength=n_bins)
    c = np.bincount(idx[ok], weights=(mass * np.cos(2 * ph))[ok], minlength=n_bins)
    s = np.bincount(idx[ok], weights=(mass * np.sin(2 * ph))[ok], minlength=n_bins)
    return 0.5 * (e[1:] + e[:-1]), np.hypot(c, s) / np.maximum(a0, 1e-300)


def clump_fraction(snap, mass, pixel=CLUMP_PIXEL, contrast=CLUMP_CONTRAST, R_min=CLUMP_R_MIN,
                   half_width=MAP_HALF_WIDTH):
    """Topaklanma: R_min < R < half_width'te, aynı yarıçaptaki azimut ortalamasının
    `contrast` katından yoğun piksellerdeki kütle kesri. Parçalanmanın ölçüsü. Boyutsuz.
    """
    L = half_width
    nb = int(round(2 * L / pixel))
    H = np.histogram2d(snap[:, 0], snap[:, 1], nb, [[-L, L], [-L, L]], weights=mass)[0]
    c = -L + pixel * (np.arange(nb) + 0.5)
    X, Y = np.meshgrid(c, c, indexing="ij")
    Rp = np.hypot(X, Y)
    ring = np.floor(Rp / (2 * pixel)).astype(int)
    ring_mean = np.bincount(ring.ravel(), weights=H.ravel()) / np.bincount(ring.ravel())
    mean = ring_mean[ring]
    ok = (Rp > R_min) & (Rp < L) & (mean > 0)
    return float(H[ok & (H > contrast * mean)].sum() / H[ok].sum())


def summary_table(run: Run, times=(0.0, 500.0, 1000.0, 2000.0)):
    """Seçili zamanlarda R_d, kalınlık, A2/A0, topak kesri. dönüş : list[dict]."""
    rows = []
    for T in times:
        if T > run.t_snap[-1] + 1e-6:
            continue
        k = run.index_at(T)
        s = run.snapshot(k)
        rows.append(dict(t=float(run.t_snap[k]), R_d=scale_length(s, run.mass),
                         thickness=disk_thickness(s), A2=bar_strength_snapshot(s, run.mass),
                         clumps=clump_fraction(s, run.mass)))
    return rows


# ==========================================================================
# Çubuk: A2/A0(R), m=2 fazı, çubuk uzunluğu
# ==========================================================================

BAR_DR = 0.5               # kpc, radyal kuşak genişliği
BAR_PROFILE_R_MAX = 10.0   # kpc, m=2 profilinin dış sınırı
BAR_MIN_COUNT = 200        # bir kuşağın max/faz hesabına girmesi için en az parçacık
BAR_PHASE_TOL = np.deg2rad(10.0)  # faz sabitliği toleransı (çubuk uzunluğu)
INNER_R = 3.0              # kpc, "merkezi bölge" kalınlığı için (bükülme testi)
BAR_MIN_LENGTH = 2.0       # kpc, faz-sabit bölge en az bu kadar uzunsa m=2 yapısı "çubuk" sayılır
BAR_PERSIST = 100.0        # Myr, çubuk oluşum zamanı: çubuk ölçütü bu süre kesintisiz sağlanmalı


def m2_profile(snap, mass, dr=BAR_DR, R_max=BAR_PROFILE_R_MAX):
    """Radyal kuşaklarda m=2 Fourier bileşeni.

    dönüş : dict — R_mid [kpc], A2 (A2/A0, boyutsuz), phase [rad, (-pi/2, pi/2]], count
    phase = (1/2) atan2(sum m sin 2phi, sum m cos 2phi): çubuğun uzun ekseninin açısı.
    """
    x, y = snap[:, 0], snap[:, 1]
    R, ph = np.hypot(x, y), np.arctan2(y, x)
    e = np.arange(0.0, R_max + 0.5 * dr, dr)
    nb = len(e) - 1
    idx = np.digitize(R, e) - 1
    ok = (idx >= 0) & (idx < nb)
    a0 = np.bincount(idx[ok], weights=mass[ok], minlength=nb)
    c = np.bincount(idx[ok], weights=(mass * np.cos(2 * ph))[ok], minlength=nb)
    s = np.bincount(idx[ok], weights=(mass * np.sin(2 * ph))[ok], minlength=nb)
    return dict(R_mid=0.5 * (e[1:] + e[:-1]), A2=np.hypot(c, s) / np.maximum(a0, 1e-300),
                phase=0.5 * np.arctan2(s, c), count=np.bincount(idx[ok], minlength=nb), edges=e)


def _dphase(a, b):
    """m=2 faz farkı, (-pi/2, pi/2] aralığına sarılmış."""
    return (a - b + 0.5 * np.pi) % np.pi - 0.5 * np.pi


def bar_measure(snap, mass, threshold=BAR_THRESHOLD, tol=BAR_PHASE_TOL, min_count=BAR_MIN_COUNT,
                R_search=A2_R_MAX):
    """m=2 gücü, çubuk gücü ve çubuk uzunluğu.

    A2max_all, R_peak_all : tüm kuşaklardaki A2/A0(R) maksimumu ve yarıçapı — sarmal kolları da yakalar
    A2max, R_peak         : R < R_search kuşaklarındaki maksimum (çubuk gücü); bölge, toplu A2'nin
                            ölçüldüğü bölgeyle aynıdır, yalnızca radyal olarak çözülmüştür
    L_bar : tepe kuşağından dışarı doğru m=2 fazının tepe fazından en fazla `tol` saptığı son
            kuşağın dış kenarı [kpc]; A2max < threshold ise NaN (çubuk yok)
    """
    p = m2_profile(snap, mass)
    good = p["count"] >= min_count
    A2_all = np.where(good, p["A2"], -1.0)
    i_all = int(np.argmax(A2_all))
    A2_in = np.where(good & (p["edges"][1:] <= R_search + 1e-9), p["A2"], -1.0)
    i0 = int(np.argmax(A2_in))
    A2max = float(p["A2"][i0])
    L = np.nan
    if A2max >= threshold:
        j = i0
        while j + 1 < len(A2_all) and good[j + 1] and abs(_dphase(p["phase"][j + 1], p["phase"][i0])) < tol:
            j += 1
        L = float(p["edges"][j + 1])
    return dict(A2max=A2max, R_peak=float(p["R_mid"][i0]), L_bar=L, phase=float(p["phase"][i0]),
                A2max_all=float(p["A2"][i_all]), R_peak_all=float(p["R_mid"][i_all]))


def bar_present(A2max, L_bar, threshold=BAR_THRESHOLD, min_length=BAR_MIN_LENGTH):
    """Çubuk ölçütü: A2max >= threshold ve faz-sabit uzunluk >= min_length. Dizi ya da skaler."""
    L = np.nan_to_num(np.asarray(L_bar, dtype=float), nan=0.0)
    return (np.asarray(A2max) >= threshold) & (L >= min_length)


def bar_onset_time(t, present, persist=BAR_PERSIST):
    """Çubuk ölçütünün `persist` Myr kesintisiz sağlandığı ilk aralığın başlangıcı [Myr]; yoksa None."""
    t = np.asarray(t)
    start = None
    for ti, ok in zip(t, present):
        if ok:
            start = ti if start is None else start
            if ti - start >= persist - 1e-9:
                return float(start)
        else:
            start = None
    return None


def inner_thickness(snap, R_max=INNER_R):
    """Merkezi bölgenin (R < R_max) dikey kalınlığı: medyan |z - <z>| [kpc]."""
    R = np.hypot(snap[:, 0], snap[:, 1])
    z = snap[R < R_max, 2]
    return float(np.median(np.abs(z - z.mean())))


def bar_history(run: Run, use_cache=True):
    """Her snapshot için bar_measure çıktıları, iç kalınlık ve toplu A2 (R<A2_R_MAX).

    Sonuç runs/<ad>/bar_history.npz'ye önbelleğe alınır (snapshot'lar değişmez).
    dönüş : dict, diziler (S,)
    """
    cache = run.path / "bar_history.npz"
    keys = ("t", "A2max", "R_peak", "L_bar", "phase", "A2max_all", "R_peak_all", "z_inner", "A2_global")
    if use_cache and cache.exists():
        d = dict(np.load(cache))
        if all(k in d for k in keys) and len(d["t"]) == run.n_snap:
            return d
    out = {k: np.empty(run.n_snap) for k in keys}
    for k in range(run.n_snap):
        s = run.snapshot(k)
        b = bar_measure(s, run.mass)
        out["t"][k] = run.t_snap[k]
        for key in ("A2max", "R_peak", "L_bar", "phase", "A2max_all", "R_peak_all"):
            out[key][k] = b[key]
        out["z_inner"][k] = inner_thickness(s)
        out["A2_global"][k] = bar_strength_snapshot(s, run.mass)
    np.savez(cache, **out)
    return out


# ==========================================================================
# Desen hızı, eşdönme, merkezi kütle
# ==========================================================================

PATTERN_FIT_HALF = 5        # snapshot; Omega_p için ±5 snapshot (±50 Myr) doğrusal uyum
PATTERN_T_MIN = 500.0       # Myr; bundan önce çubuk oluşum evresi, faz açması güvenilmez
PATTERN_R_IN = 0.5          # kpc; çubuk fazı R_IN < R < L_bar bölgesinden ölçülür
CR_RING_N = 120             # eşdönme için v_c ölçüm halkası sayısı
CR_RING_R = (0.25, 15.0)    # kpc
CR_VC_SMOOTH = 0.5          # kpc
OBS_CR_RATIO = (1.0, 1.4)   # gözlenen çubuklarda R_CR / R_çubuk aralığı
CENTRAL_R = 1.0             # kpc, merkezi kütle yarıçapı


def bar_region_phase(snap, mass, bm):
    """Çubuk bölgesinin (R_IN < R < L_bar; çubuk yoksa R_peak + kuşak) m=2 fazı [rad]."""
    x, y = snap[:, 0], snap[:, 1]
    R, ph = np.hypot(x, y), np.arctan2(y, x)
    R_out = bm["L_bar"] if np.isfinite(bm["L_bar"]) else bm["R_peak"] + BAR_DR
    sel = (R > PATTERN_R_IN) & (R < R_out)
    return 0.5 * np.arctan2(np.sum(mass[sel] * np.sin(2 * ph[sel])), np.sum(mass[sel] * np.cos(2 * ph[sel])))


def half_max_length(snap, mass, bm):
    """Alternatif çubuk uzunluğu: A2/A0(R)'nin tepeden dışarı doğru tepe değerinin yarısına düştüğü yarıçap [kpc]."""
    p = m2_profile(snap, mass)
    j = int(np.argmin(np.abs(p["R_mid"] - bm["R_peak"])))
    while j + 1 < len(p["A2"]) and p["A2"][j + 1] >= 0.5 * bm["A2max"]:
        j += 1
    return float(p["edges"][j + 1])


def circular_velocity_snapshot(run: Run, k, R):
    """k. snapshot'ın gerçek alanında (parçacıklar + sabit bileşenler) v_c(R) [kpc/Myr], R : kpc."""
    from .ics import _gaussian_smooth, particle_vc2
    gp, _, num = params_from_config(run)
    s = run.snapshot(k)
    pos = np.ascontiguousarray(s[:, :3] - np.average(s[:, :3], axis=0, weights=run.mass))
    Rr = np.linspace(*CR_RING_R, CR_RING_N)
    vc2 = _gaussian_smooth(Rr, particle_vc2(gp, num, pos, run.mass, Rr), CR_VC_SMOOTH)
    return np.sqrt(np.interp(R, Rr, vc2))


def corotation_radius(R, v_c, omega_p):
    """v_c(R_CR) = Omega_p R_CR çözümü (ilk kesişim). R : kpc, v_c : kpc/Myr, omega_p : 1/Myr -> kpc ya da NaN."""
    d = v_c / R - omega_p
    i = np.flatnonzero((d[:-1] > 0) & (d[1:] <= 0))
    if not i.size:
        return np.nan
    i = i[0]
    return float(R[i] + d[i] / (d[i] - d[i + 1]) * (R[i + 1] - R[i]))


def pattern_history(run: Run, use_cache=True):
    """Her snapshot: çubuk fazı (açılmış), Omega_p [1/Myr], R_CR [kpc], L_faz, L_yarı [kpc], M_merkez [M_sun].

    Omega_p, açılmış fazın ±PATTERN_FIT_HALF snapshot penceresindeki doğrusal eğimi;
    t < PATTERN_T_MIN için NaN. m=2 fazının periyodu pi olduğundan açma, adım başına
    < 90° dönmeyi varsayar (t >= PATTERN_T_MIN'de sağlandığı kontrol edildi).
    """
    cache = run.path / "pattern_history.npz"
    keys = ("t", "phase", "omega_p", "R_CR", "L_phase", "L_half", "present", "M_central", "max_step_deg")
    if use_cache and cache.exists():
        d = dict(np.load(cache))
        if all(k in d for k in keys) and len(d["t"]) == run.n_snap:
            return d
    n = run.n_snap
    t = run.t_snap
    phase = np.empty(n)
    L_half = np.empty(n)
    M_c = np.empty(n)
    bms = []
    for k in range(n):
        s = run.snapshot(k)
        bm = bar_measure(s, run.mass)
        bms.append(bm)
        phase[k] = bar_region_phase(s, run.mass, bm)
        L_half[k] = half_max_length(s, run.mass, bm)
        M_c[k] = run.mass[np.hypot(s[:, 0], s[:, 1]) < CENTRAL_R].sum()
    un = np.unwrap(2.0 * phase) / 2.0
    steps = np.degrees(np.abs(np.diff(un)))
    W = PATTERN_FIT_HALF
    omega = np.full(n, np.nan)
    R_cr = np.full(n, np.nan)
    Rg = np.linspace(*CR_RING_R, CR_RING_N)
    for k in range(W, n - W):
        if t[k] < PATTERN_T_MIN:
            continue
        omega[k] = np.polyfit(t[k - W:k + W + 1], un[k - W:k + W + 1], 1)[0]
        R_cr[k] = corotation_radius(Rg, circular_velocity_snapshot(run, k, Rg), omega[k])
    L_phase = np.array([b["L_bar"] for b in bms])
    present = bar_present(np.array([b["A2max"] for b in bms]), L_phase)
    late = t >= PATTERN_T_MIN
    out = dict(t=t, phase=un, omega_p=omega, R_CR=R_cr, L_phase=L_phase, L_half=L_half,
               present=present, M_central=M_c,
               max_step_deg=np.array([steps[late[1:]].max() if np.any(late[1:]) else np.nan]))
    np.savez(cache, **out)
    return out


def corotation_stats(ph):
    """Çubuk varken (t >= PATTERN_T_MIN) R_CR/L oranları: iki uzunluk tanımı için medyan, çeyrekler, gözlem aralığındaki kesir."""
    sel = ph["present"] & np.isfinite(ph["R_CR"])
    out = {}
    for key in ("L_phase", "L_half"):
        r = ph["R_CR"][sel] / ph[key][sel]
        out[key] = dict(n=int(r.size), median=float(np.median(r)) if r.size else None,
                        q25=float(np.percentile(r, 25)) if r.size else None,
                        q75=float(np.percentile(r, 75)) if r.size else None,
                        frac_obs=float(np.mean((r >= OBS_CR_RATIO[0]) & (r <= OBS_CR_RATIO[1]))) if r.size else None)
    return out


# ==========================================================================
# Toomre Q(R) — parçacıkların gerçek alanından
# ==========================================================================

Q_EDGES = np.arange(1.0, 12.01, 1.0)   # kpc, Q(R) halka sınırları
Q_RING_N = 60                          # v_c ölçüm halkası sayısı
Q_VC_SMOOTH = 0.75                     # kpc, v_c^2 Gauss düzgünleştirme genişliği


def params_from_config(run: Run):
    """config.json'dan (GalaxyParams, DiskICParams, NumericalParams)."""
    from .config import (BlackHoleParams, BulgeParams, DiskICParams, DiskParams, GalaxyParams,
                         HaloParams, NumericalParams)
    c = run.config
    gp = GalaxyParams(disk=DiskParams(**c["galaxy"]["disk"]), bulge=BulgeParams(**c["galaxy"]["bulge"]),
                      halo=HaloParams(**c["galaxy"]["halo"]), bh=BlackHoleParams(**c["galaxy"]["bh"]))
    return gp, DiskICParams(**c["disk_ic"]), NumericalParams(**c["numerical"])


def toomre_Q_profile(run: Run, k, edges=Q_EDGES):
    """k. snapshot'ta Q(R) = sigma_R kappa / (3.36 G Sigma).

    Sigma ve sigma_R parçacıklardan; kappa, parçacıklar + sabit bileşenlerin halkalar üzerinde
    ölçülen v_c'sinden (ics.particle_vc2, aynı eps): kappa^2 = R dOmega^2/dR + 4 Omega^2.
    Çubuk bölgesinde sigma_R akış hareketlerini de içerir; orada Q yalnızca kabaca anlamlıdır.
    dönüş : dict — R_mid [kpc], Q, sigma_R [kpc/Myr], Sigma [M_sun/kpc^2], kappa [1/Myr], v_c [kpc/Myr]
    """
    from .ics import TOOMRE_STELLAR, DiskParticles, measure_profiles, particle_vc2, _gaussian_smooth
    from .units import G

    gp, _, num = params_from_config(run)
    s = run.snapshot(k)
    pos = np.ascontiguousarray(s[:, :3])
    pos = pos - np.average(pos, axis=0, weights=run.mass)  # kütle merkezine göre
    vel = s[:, 3:] - np.average(s[:, 3:], axis=0, weights=run.mass)
    m = measure_profiles(DiskParticles(pos, vel, run.mass), edges)

    R_ring = np.linspace(0.5 * edges[0], edges[-1] + 1.0, Q_RING_N)
    vc2 = _gaussian_smooth(R_ring, particle_vc2(gp, num, pos, run.mass, R_ring), Q_VC_SMOOTH)
    om2 = vc2 / R_ring ** 2
    k2 = R_ring * np.gradient(om2, R_ring) + 4.0 * om2
    kap = np.sqrt(np.maximum(np.interp(m["R_mid"], R_ring, k2), 0.0))
    Q = m["sigma_R"] * kap / (TOOMRE_STELLAR * G * m["Sigma"])
    return dict(R_mid=m["R_mid"], Q=Q, sigma_R=m["sigma_R"], Sigma=m["Sigma"], kappa=kap,
                v_c=np.sqrt(np.interp(m["R_mid"], R_ring, vc2)))


# ==========================================================================
# Görselleştirme (koşu bittikten sonra)
# ==========================================================================

def face_on_map(ax, snap, mass, half_width=MAP_HALF_WIDTH, bins=300, vmin=None, vmax=None):
    """Yüzey yoğunluğu haritası (log). ax : matplotlib ekseni."""
    from matplotlib.colors import LogNorm
    L = half_width
    H, xe, ye = np.histogram2d(snap[:, 0], snap[:, 1], bins, [[-L, L], [-L, L]], weights=mass)
    area = (2 * L / bins) ** 2
    S = np.where(H > 0, H / area, np.nan)
    im = ax.imshow(S.T, origin="lower", extent=[-L, L, -L, L], cmap="magma",
                   norm=LogNorm(vmin=vmin, vmax=vmax))
    ax.set(xlabel="x [kpc]", ylabel="y [kpc]", aspect="equal")
    return im


def map_color_limits(snap, mass, half_width=MAP_HALF_WIDTH, bins=300):
    """Animasyon/karşılaştırma için sabit renk aralığı (t=0 haritasından), M_sun/kpc^2."""
    L = half_width
    H = np.histogram2d(snap[:, 0], snap[:, 1], bins, [[-L, L], [-L, L]], weights=mass)[0]
    S = H[H > 0] / (2 * L / bins) ** 2
    return float(np.percentile(S, 5)), float(S.max())


def comparison_figure(run: Run, out_path, times=(0.0, 500.0, 1000.0, 2000.0), dpi=110):
    """Dört kareli (yüzden) karşılaştırma PNG'si."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    vmin, vmax = map_color_limits(run.snapshot(0), run.mass)
    fig, axes = plt.subplots(1, len(times), figsize=(4.2 * len(times), 4.6))
    for ax, T in zip(axes, times):
        k = run.index_at(T)
        s = run.snapshot(k)
        face_on_map(ax, s, run.mass, vmin=vmin, vmax=vmax)
        b = bar_measure(s, run.mass)
        L = f"{b['L_bar']:.1f} kpc" if bar_present(b["A2max"], b["L_bar"]) else "çubuk yok"
        ax.set_title(f"t = {run.t_snap[k]:.0f} Myr   A2max = {b['A2max']:.2f}   L = {L}", fontsize=10)
    fig.suptitle(run_label(run))
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def run_label(run: Run):
    c = run.config
    return (f"{run.name}: N={c['numerical']['N']}, Q={c['disk_ic']['Q_target']}, "
            f"M_vir={c['galaxy']['halo']['M_vir']:.2g} M☉, ELN={c['eln']:.2f}, "
            f"öz-yerçekimi={c['self_gravity']}")


def make_animation(run: Run, out_path, fps=12, every=1, dpi=80):
    """Snapshot'lardan yüzden yoğunluk + A2/A0(t) animasyonu (MP4: ffmpeg, GIF: Pillow)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FFMpegWriter, PillowWriter

    vmin, vmax = map_color_limits(run.snapshot(0), run.mass)
    h = bar_history(run)
    t_ts, A2 = h["t"], h["A2max"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw=dict(width_ratios=[1, 1.15]))
    bx.plot(t_ts, A2, color="0.75", lw=1)
    (trace,) = bx.plot([], [], color="#d95f02", lw=1.6)
    marker = bx.axvline(0.0, color="k", lw=0.8)
    bx.axhline(BAR_THRESHOLD, color="0.5", ls=":", lw=1)
    bx.set(xlabel="t [Myr]", ylabel="maks A2/A0(R)",
           title=f"Çubuk gücü: A2/A0(R) maksimumu, R < {A2_R_MAX:g} kpc",
           xlim=(0, t_ts[-1]), ylim=(0, max(0.05, 1.1 * A2.max())))
    bx.grid(alpha=0.25)
    fig.suptitle(run_label(run), fontsize=9)

    writer = FFMpegWriter(fps=fps, bitrate=2400) if str(out_path).endswith(".mp4") else PillowWriter(fps=fps)
    with writer.saving(fig, str(out_path), dpi):
        for k in range(0, run.n_snap, every):
            ax.clear()
            face_on_map(ax, run.snapshot(k), run.mass, vmin=vmin, vmax=vmax)
            t = run.t_snap[k]
            ax.set_title(f"t = {t:.0f} Myr")
            sel = t_ts <= t + 1e-9
            trace.set_data(t_ts[sel], A2[sel])
            marker.set_xdata([t, t])
            writer.grab_frame()
    plt.close(fig)


def bar_overview_figure(runs, out_path, dpi=110):
    """Üst: çubuk gücü (A2/A0(R) maksimumu, R < A2_R_MAX) ve ince çizgiyle tüm-R m=2 maksimumu.
    Alt: çubuk uzunluğu (faz-sabit bölge), yalnızca çubuk ölçütü sağlanırken."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    for i, r in enumerate(runs):
        h = bar_history(r)
        c = f"C{i}"
        ls = "--" if r.name.startswith("dogrulama") else "-"
        lab = f"{r.name} (N={r.config['numerical']['N']}, {r.config['self_gravity']})"
        ax.plot(h["t"], running_mean(h["t"], h["A2max"], 30.0), ls, color=c, lw=1.4, label=lab)
        ax.plot(h["t"], h["A2max_all"], ls, color=c, lw=0.5, alpha=0.5)
        pres = bar_present(h["A2max"], h["L_bar"])
        bx.plot(h["t"], np.where(pres, h["L_bar"], np.nan), "o" if ls == "-" else "x", color=c, ms=2.5, label=r.name)
    ax.axhline(BAR_THRESHOLD, color="0.5", ls=":", lw=1)
    ax.set(ylabel="maks A2/A0(R)", title=f"Çubuk gücü: kalın — R < {A2_R_MAX:g} kpc kuşaklarında maksimum (30 Myr ort.); "
           "ince — tüm yarıçaplarda maksimum (sarmal kollar dahil)")
    bx.set(xlabel="t [Myr]", ylabel="L_bar [kpc]",
           title=f"Çubuk uzunluğu: m=2 fazının ±{np.degrees(BAR_PHASE_TOL):.0f}° sabit kaldığı bölge "
                 f"(yalnızca A2max ≥ {BAR_THRESHOLD} ve L ≥ {BAR_MIN_LENGTH:g} kpc iken)")
    for a in (ax, bx):
        a.grid(alpha=0.25)
    ax.legend(fontsize=8)
    ax.title.set_fontsize(9)
    bx.title.set_fontsize(9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def m2_profiles_figure(runs, out_path, dpi=100):
    """Her koşu için seçili zamanlarda A2/A0(R) ve m=2 fazı(R). Faz sabit bölge = çubuk."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(runs), 2, figsize=(12, 3.0 * len(runs)), squeeze=False)
    for (a, b), r in zip(axes, runs):
        for j, T in enumerate(panel_times(r)[1:]):
            k = r.index_at(T)
            s = r.snapshot(k)
            p = m2_profile(s, r.mass)
            bm = bar_measure(s, r.mass)
            ok = p["count"] >= BAR_MIN_COUNT
            lab = f"t={r.t_snap[k]:.0f}"
            if bar_present(bm["A2max"], bm["L_bar"]):
                lab += f", L={bm['L_bar']:.1f} kpc"
                b.axvspan(0, bm["L_bar"], color=f"C{j}", alpha=0.07)
            a.plot(p["R_mid"][ok], p["A2"][ok], "-o", ms=2.5, color=f"C{j}", label=lab)
            b.plot(p["R_mid"][ok], np.degrees(_dphase(p["phase"][ok], bm["phase"])), "-o", ms=2.5, color=f"C{j}")
        a.axhline(BAR_THRESHOLD, color="0.5", ls=":", lw=1)
        a.axvline(A2_R_MAX, color="0.5", ls="--", lw=0.8)
        a.set(ylabel="A2/A0(R)", title=f"{r.name}: A2/A0(R)")
        b.axhspan(-np.degrees(BAR_PHASE_TOL), np.degrees(BAR_PHASE_TOL), color="0.85")
        b.set(ylabel="faz − çubuk fazı [°]", ylim=(-95, 95), title=f"{r.name}: m=2 fazı (gri bant ±{np.degrees(BAR_PHASE_TOL):.0f}°)")
        a.legend(fontsize=7)
        for x in (a, b):
            x.grid(alpha=0.25)
            x.set_xlabel("R [kpc]")
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


Q_TRACK_TIMES = (0.0, 200.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1200.0, 1500.0, 2000.0)
Q_TRACK_BANDS = ((2.0, 4.0), (4.0, 8.0), (8.0, 12.0))   # kpc


def q_history(run: Run, times=Q_TRACK_TIMES, bands=Q_TRACK_BANDS):
    """Seçili zamanlarda halka bantlarında ortalama Q. dönüş : (t [Myr], Q (len(times), len(bands)))."""
    ts, Qs = [], []
    for T in times:
        if T > run.t_snap[-1] + 1e-6:
            continue
        k = run.index_at(T)
        q = toomre_Q_profile(run, k)
        ts.append(float(run.t_snap[k]))
        Qs.append([np.mean(q["Q"][(q["R_mid"] > lo) & (q["R_mid"] < hi)]) for lo, hi in bands])
    return np.array(ts), np.array(Qs)


def bar_weakening_figure(runs, out_path, dpi=100):
    """Hipotez testi. Sol: A2max(t) ve iç kalınlık z(R<INNER_R)(t) (bükülme).
    Sağ: A2max(t) ve Q(t) bantları (ısınma)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(runs), 2, figsize=(13, 3.6 * len(runs)), squeeze=False)
    for (a, b), r in zip(axes, runs):
        h = bar_history(r)
        A = running_mean(h["t"], h["A2max"], 30.0)
        for ax in (a, b):
            ax.plot(h["t"], A, color="k", lw=1.4, label="çubuk gücü A2max (30 Myr ort.)")
            ax.set(xlabel="t [Myr]", ylabel="A2max")
            ax.grid(alpha=0.25)
        a2 = a.twinx()
        a2.plot(h["t"], h["z_inner"] * 1e3, color="C3", lw=1.2, label=f"medyan|z|, R < {INNER_R:g} kpc")
        a2.set_ylabel("iç kalınlık [pc]", color="C3")
        a.set_title(f"{r.name}: bükülme testi — A2 düşerken iç kalınlık sıçrıyor mu?", fontsize=9)
        tq, Q = q_history(r)
        b2 = b.twinx()
        for j, (lo, hi) in enumerate(Q_TRACK_BANDS):
            b2.plot(tq, Q[:, j], "-o", ms=3, color=f"C{j}", label=f"Q, {lo:g}–{hi:g} kpc")
        b2.set_ylabel("Toomre Q")
        b.set_title(f"{r.name}: ısınma testi — çubuk zayıflarken Q yükseliyor mu?", fontsize=9)
        h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
        a.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper left")
        h1, l1 = b.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels()
        b.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def q_profile_figure(run: Run, out_path, times=(0.0, 200.0, 800.0), target=None, dpi=110):
    """Q(R) profillerinin karşılaştırması (öz-düzenleme)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for T in times:
        k = run.index_at(T)
        q = toomre_Q_profile(run, k)
        ax.plot(q["R_mid"], q["Q"], "-o", ms=4, label=f"t = {run.t_snap[k]:.0f} Myr (ort. {np.mean(q['Q']):.2f})")
    if target is not None:
        ax.axhline(target, color="k", ls="--", lw=1, label=f"Q_target = {target}")
    ax.axhline(1.0, color="crimson", ls=":", lw=1, label="Q = 1")
    ax.set(xlabel="R [kpc]", ylabel="Toomre Q", title=f"{run.name}: Q(R) — öz-düzenleme", ylim=(0, None))
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def pattern_speed_figure(runs, out_path, window=(1300.0, 1800.0), dpi=100):
    """Üst: Omega_p(t). Alt: R_CR / L_çubuk (iki uzunluk tanımı), gözlenen 1.0–1.4 bandı."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    for i, r in enumerate(runs):
        ph = pattern_history(r)
        c = f"C{i}"
        ax.plot(ph["t"], ph["omega_p"] / KMS_KPMY, color=c, lw=1.4, label=r.name)
        sel = ph["present"] & np.isfinite(ph["R_CR"])
        bx.plot(ph["t"][sel], (ph["R_CR"] / ph["L_phase"])[sel], "o", ms=3, color=c, label=f"{r.name}: L = faz-sabit")
        bx.plot(ph["t"][sel], (ph["R_CR"] / ph["L_half"])[sel], "x", ms=3, color=c, alpha=0.6, label=f"{r.name}: L = A2 yarı-maks")
    ax.axvspan(*window, color="0.9", zorder=0)
    ax.set(ylabel="Ω_p [km/s/kpc]", title="Çubuk desen hızı (m=2 fazının ±50 Myr doğrusal eğimi); gri: ana koşunun zayıflama dönemi")
    bx.axhspan(*OBS_CR_RATIO, color="C2", alpha=0.12, label=f"gözlenen {OBS_CR_RATIO[0]}–{OBS_CR_RATIO[1]}")
    bx.set(xlabel="t [Myr]", ylabel="R_CR / L_çubuk", ylim=(0, 3.5),
           title="Eşdönme / çubuk uzunluğu (yalnızca çubuk ölçütü sağlanırken)")
    for a in (ax, bx):
        a.grid(alpha=0.25)
        a.title.set_fontsize(9)
    ax.legend(fontsize=8)
    bx.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def central_mass_figure(runs, out_path, window=(1300.0, 1800.0), dpi=100):
    """Her koşu: çubuk gücü A2max(t) ve R < CENTRAL_R içindeki disk kütlesi (aynı grafikte, ikiz eksen)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(runs), 1, figsize=(11, 3.4 * len(runs)), squeeze=False)
    for ax, r in zip(axes[:, 0], runs):
        h = bar_history(r)
        ph = pattern_history(r)
        ax.plot(h["t"], running_mean(h["t"], h["A2max"], 30.0), color="k", lw=1.4, label="çubuk gücü A2max (30 Myr ort.)")
        ax.set(ylabel="A2max", xlabel="t [Myr]")
        ax.axvspan(*window, color="0.92", zorder=0)
        a2 = ax.twinx()
        a2.plot(ph["t"], ph["M_central"] / 1e9, color="C3", lw=1.4, label=f"disk kütlesi, R < {CENTRAL_R:g} kpc")
        a2.set_ylabel("M [10⁹ M☉]", color="C3")
        ax.set_title(f"{r.name}: merkezi kütle yoğunlaşması ve çubuk gücü", fontsize=9)
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper left")
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


# ==========================================================================
# Toplu çıktı üretimi:  python -m src.analysis
# ==========================================================================

MAIN_RUNS = ("ana", "soguk", "hafif_hale")
CHECK_RUNS = ("dogrulama_dogrudan", "dogrulama_agac", "dogrulama_dogrudan_seed7")
ANIMATE = ("ana", "hafif_hale")
PANEL_TIMES = (0.0, 500.0, 1000.0, 2000.0)


def panel_times(run: Run):
    """İstenen 0/500/1000/2000 Myr; daha kısa koşularda t_end'e eşit aralıklı dört an."""
    if run.t_snap[-1] >= PANEL_TIMES[-1] - 1e-6:
        return PANEL_TIMES
    return tuple(np.linspace(0.0, run.t_snap[-1], len(PANEL_TIMES)))


def a2_overlay_figure(runs, out_path, dpi=110):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(11, 5))
    for r in runs:
        t, A2 = bar_series(r)
        ls = "--" if r.name.startswith("dogrulama") else "-"
        ax.plot(t, A2, ls, lw=1.3 if ls == "-" else 1.0,
                label=f"{r.name} (N={r.config['numerical']['N']}, {r.config['self_gravity']})")
    ax.axhline(BAR_THRESHOLD, color="0.5", ls=":", lw=1)
    ax.set(xlabel="t [Myr]", ylabel="A2/A0", title=f"Çubuk gücü A2/A0 (R < {A2_R_MAX:g} kpc)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def run_metrics(run: Run):
    t, A2 = bar_series(run)
    _, dE = energy_error(run)
    late = t >= 0.75 * t[-1]
    r_com, v_com = com_drift(run)
    h = bar_history(run)
    pres = bar_present(h["A2max"], h["L_bar"])
    hl = h["t"] >= 0.75 * h["t"][-1]
    Lp = np.where(pres, h["L_bar"], np.nan)
    bar = dict(t_onset=bar_onset_time(h["t"], pres), present_frac=float(pres.mean()),
               present_frac_late=float(pres[hl].mean()), A2max_peak=float(h["A2max"].max()),
               t_A2max_peak=float(h["t"][h["A2max"].argmax()]), A2max_late_mean=float(h["A2max"][hl].mean()),
               L_bar_median=float(np.nanmedian(Lp)) if np.any(pres) else None,
               L_bar_late_median=float(np.nanmedian(Lp[hl])) if np.any(pres[hl]) else None,
               A2max_t0=float(h["A2max"][0]))
    if run.name in MAIN_RUNS:
        ph = pattern_history(run)
        bar["corotation"] = corotation_stats(ph)
        bar["omega_p_kms_kpc"] = {str(int(T)): float(ph["omega_p"][run.index_at(T)] / KMS_KPMY)
                                  for T in (500, 1000, 1300, 1500, 1800) if T <= run.t_snap[-1]}
        bar["max_phase_step_deg_after_tmin"] = float(ph["max_step_deg"][0])
        bar["M_central"] = {str(int(T)): float(ph["M_central"][run.index_at(T)]) for T in (0, 500, 1300, 1800, 2000)}
    return dict(bar=bar,
        name=run.name, N=run.config["numerical"]["N"], Q=run.config["disk_ic"]["Q_target"],
        M_vir=run.config["galaxy"]["halo"]["M_vir"], eln=run.config["eln"],
        self_gravity=run.config["self_gravity"], t_end=float(t[-1]),
        t_bar=bar_formation_time(t, A2), A2_max=float(A2.max()), t_A2_max=float(t[A2.argmax()]),
        A2_late_mean=float(A2[late].mean()), A2_late_std=float(A2[late].std()),
        dE_max=float(np.abs(dE).max()), dE_end=float(dE[-1]), dE_rate_per_Myr=energy_drift_rate(run),
        com_max_kpc=r_com, vcom_max_kms=v_com / KMS_KPMY,
        table=summary_table(run, panel_times(run)),
    )


def main(runs_dir="runs", animate=True):
    root = Path(runs_dir)
    present = [n for n in MAIN_RUNS + CHECK_RUNS if is_complete(root / n)]
    missing = [n for n in MAIN_RUNS + CHECK_RUNS if n not in present]
    if missing:
        print(f"[analiz] eksik/bitmemiş koşular atlandı: {missing}")
    runs = [load_run(n, root) for n in present]
    a2_overlay_figure(runs, root / "a2_toplu_R5.png")   # eski tanım (R<5 toplu), karşılaştırma için
    bar_overview_figure(runs, root / "a2_karsilastirma.png")
    profile_runs = [r for r in runs if r.name in MAIN_RUNS + ("dogrulama_dogrudan", "dogrulama_agac")]
    m2_profiles_figure(profile_runs, root / "cubuk_faz_profilleri.png")
    weak_runs = [r for r in runs if r.name in ("ana", "hafif_hale")]
    if weak_runs:
        bar_weakening_figure(weak_runs, root / "cubuk_zayiflama_testi.png")
    main_runs = [r for r in runs if r.name in MAIN_RUNS]
    if main_runs:
        pattern_speed_figure(main_runs, root / "desen_hizi.png")
        central_mass_figure(main_runs, root / "merkezi_kutle.png")
    if "soguk" in present:
        r = next(r for r in runs if r.name == "soguk")
        q_profile_figure(r, r.path / "Q_profili.png", target=r.config["disk_ic"]["Q_target"])
    metrics = []
    for r in runs:
        comparison_figure(r, r.path / "karsilastirma.png", times=panel_times(r))
        metrics.append(run_metrics(r))
        b = metrics[-1]["bar"]
        print(f"[analiz] {r.name}: çubuk oluşumu {b['t_onset']} Myr, maks çubuk gücü {b['A2max_peak']:.3f}, "
              f"çubuk uzunluğu (medyan) {b['L_bar_median']} kpc, maks|dE/E| {metrics[-1]['dE_max']:.1e}", flush=True)
    (root / "ozet_metrikler.json").write_text(json.dumps(metrics, indent=2, ensure_ascii=False))
    if animate:
        for r in runs:
            if r.name in ANIMATE:
                make_animation(r, r.path / "animasyon.mp4")
                print(f"[analiz] animasyon: {r.path / 'animasyon.mp4'}", flush=True)
    return metrics


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Aşama 3 çıktıları: A2/A0 karşılaştırması, 4 kareli PNG'ler, animasyon")
    ap.add_argument("--runs_dir", default="runs")
    ap.add_argument("--no-animation", action="store_true")
    a = ap.parse_args()
    main(a.runs_dir, animate=not a.no_animation)
