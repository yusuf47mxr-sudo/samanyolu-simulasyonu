"""N-cisim disk başlangıç koşulları (Aşama 2a). Zaman evrimi yok.

Disk:
    Sigma(R) = Sigma_0 exp(-R/R_d), R < R_cut
    rho(R,z) ∝ Sigma(R) sech^2(z/z_0)

Hızlar:
    sigma_R   : Q(R) = Q_target her yarıçapta  ->  sigma_R = Q 3.36 G Sigma / kappa
    sigma_phi : epicycle yaklaşımı  sigma_phi^2 = sigma_R^2 kappa^2 / (4 Omega^2)
    sigma_z   : izotermal levha     sigma_z^2   = pi G Sigma z_0
    <v_phi>   : Jeans denklemi, asimetrik sürüklenme
        <v_phi>^2 = v_c^2 + sigma_R^2 [ dln(Sigma sigma_R^2)/dlnR + 1 - kappa^2/(4 Omega^2) ]

v_c, Omega, kappa toplam potansiyelden alınır: sabit bileşenler (kara delik,
şişkinlik, hale) + diskin kendi potansiyeli (Freeman üstel disk, analitik).

Birimler: konum : kpc, hız : kpc/Myr, kütle : M_sun, Sigma : M_sun/kpc^2.
"""

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import brentq

from .config import DiskICParams, DiskParams, GalaxyParams, NumericalParams, param_range
from .potentials import CompositePotential, FreemanDisk, build_external_potential
from .rotation import circular_velocity, kappa, omega
from .units import G, kms_to_kpcmyr, kpcmyr_to_kms

TOOMRE_STELLAR = 3.36    # yıldız diski için Toomre sabiti
ELN_BAR_THRESHOLD = 1.1  # Efstathiou–Lake–Negroponte (1982) çubuk eşiği
DLNR_STEP = 1e-4         # d ln kappa / d ln R için merkezi fark adımı (ln R cinsinden)


@dataclass
class DiskParticles:
    """pos, vel : (N,3) kpc, kpc/Myr; mass : (N,) M_sun.

    R_ring, dvc2 : v_c^2 düzeltmesi — parçacıklardan ölçülen v_c^2 eksi analitik
    (ince disk) v_c^2, düzgünleştirilmiş; kpc ve (kpc/Myr)^2.
    """
    pos: np.ndarray
    vel: np.ndarray
    mass: np.ndarray
    R_ring: np.ndarray = None
    dvc2: np.ndarray = None

    def vc2_correction(self, R):
        """R : kpc -> v_c^2 düzeltmesi : (kpc/Myr)^2."""
        return np.interp(R, self.R_ring, self.dvc2)


def ic_potential(gp: GalaxyParams, ic: DiskICParams) -> CompositePotential:
    """IC için toplam potansiyel: sabit bileşenler + üstel disk (Freeman, kütle ic.M_d)."""
    ext = build_external_potential(gp)
    return CompositePotential(ext.components + [FreemanDisk(ic.M_d, ic.R_d, ic.R_cut)])


def solve_disk_mass(gp: GalaxyParams, ic: DiskICParams, v_target_kms, R):
    """Üstel disk + şişkinlik + hale + kara delik toplamı v_c(R) = v_target verecek M_d.

    Yalnızca varsayılanı belirlemek için kullanılır (config.DiskICParams.M_d);
    arayüzde parametre değiştikçe yeniden çağrılmaz.
    v_target_kms : km/s, R : kpc  ->  M_d : M_sun
    """
    lo, hi = param_range(DiskParams, "M_d")
    v_t = kms_to_kpcmyr(v_target_kms)
    f = lambda M: float(circular_velocity(ic_potential(gp, replace(ic, M_d=M)), R)) - v_t
    return brentq(f, lo, hi, xtol=1.0, rtol=1e-12)


def eln_parameter(gp: GalaxyParams, ic: DiskICParams):
    """Efstathiou–Lake–Negroponte (1982): eps_m = V_max / sqrt(G M_d / R_d).

    V_max, IC potansiyelinin (disk + sabit bileşenler) dönme eğrisinin maksimumu.
    eps_m < 1.1 ise disk çubuk kararsızlığına açıktır. Boyutsuz.
    """
    R = np.linspace(0.05, ic.R_cut, 3000)
    v_max = circular_velocity(ic_potential(gp, ic), R).max()
    return v_max / np.sqrt(G * ic.M_d / ic.R_d)


def disk_profiles(gp: GalaxyParams, ic: DiskICParams, R):
    """Analitik radyal profiller. R : kpc (artan dizi).

    dönüş : dict — Sigma [M_sun/kpc^2], v_c, sigma_R, sigma_phi, sigma_z, vphi_mean [kpc/Myr],
            kappa, Omega [1/Myr], drift_term [boyutsuz, köşeli parantez içi]
    """
    pot = ic_potential(gp, ic)
    disk = pot["disk"]
    R = np.asarray(R, dtype=float)

    Sigma = disk.surface_density(R)
    vc = circular_velocity(pot, R)
    Om = omega(pot, R)
    ka = kappa(pot, R)

    sigR = ic.Q_target * TOOMRE_STELLAR * G * Sigma / ka
    sigphi = sigR * ka / (2.0 * Om)
    sigz = np.sqrt(np.pi * G * Sigma * ic.z_0)

    # d ln(Sigma sigma_R^2) / d ln R : sigma_R ∝ Sigma / kappa olduğundan
    #   d ln sigma_R^2 / d ln R = 2 (-R/R_d - d ln kappa / d ln R).
    # kappa'nın türevi (Phi'nin 3. türevi) her R'de merkezi farkla alınır; ızgaradan bağımsız.
    h = DLNR_STEP
    dln_kappa = (np.log(kappa(pot, R * np.exp(h))) - np.log(kappa(pot, R * np.exp(-h)))) / (2.0 * h)
    dln_sig2 = 2.0 * (-R / ic.R_d - dln_kappa)
    drift = -R / ic.R_d + dln_sig2 + 1.0 - (sigphi / sigR) ** 2
    vphi2 = vc ** 2 + sigR ** 2 * drift

    return dict(Sigma=Sigma, v_c=vc, kappa=ka, Omega=Om, sigma_R=sigR, sigma_phi=sigphi,
                sigma_z=sigz, drift_term=drift, vphi_mean=np.sqrt(np.maximum(vphi2, 0.0)),
                vphi2_negative=vphi2 < 0.0)


def sample_radii(rng, n, R_d, R_cut):
    """Kesilmiş üstel diskten R örnekle (ters CDF, ızgara üzerinde). R_d, R_cut : kpc."""
    R = np.linspace(0.0, R_cut, 20001)
    cdf = 1.0 - (1.0 + R / R_d) * np.exp(-R / R_d)
    return np.interp(rng.uniform(0.0, cdf[-1], n), cdf, R)


def sample_heights(rng, n, z_0):
    """sech^2(z/z_0) profilinden z örnekle: z = z_0 artanh(2u - 1). z_0 : kpc."""
    u = rng.uniform(0.0, 1.0, n)
    return z_0 * np.arctanh(2.0 * u - 1.0)


def particle_vc2(gp: GalaxyParams, num: NumericalParams, pos, mass, R_ring):
    """z = 0 halkalarında parçacıklar + sabit bileşenlerin verdiği v_c^2 = R <-a_R>.

    Parçacık kuvveti doğrudan toplamla, N-cisimle aynı Plummer eps ile hesaplanır.
    pos : (N,3) kpc, mass : M_sun, R_ring : kpc  ->  (kpc/Myr)^2
    """
    from .nbody import direct_accel  # nbody bu modülü içe aktarır; döngüyü önlemek için burada
    from .orbits import acceleration

    phi = np.linspace(0.0, 2.0 * np.pi, num.ic_ring_nphi, endpoint=False)
    c, s = np.tile(np.cos(phi), R_ring.size), np.tile(np.sin(phi), R_ring.size)
    Rr = np.repeat(R_ring, phi.size)
    pts = np.column_stack((Rr * c, Rr * s, np.zeros_like(Rr)))
    a = direct_accel(pts, pos, mass, num.eps, False) \
        + acceleration(build_external_potential(gp, num.eps), pts)
    aR = -(a[:, 0] * c + a[:, 1] * s).reshape(R_ring.size, phi.size).mean(axis=1)
    return R_ring * aR


def _gaussian_smooth(x, y, width):
    """y(x)'i x ekseninde Gauss çekirdeğiyle düzgünleştir. width : x biriminde."""
    w = np.exp(-0.5 * ((x[:, None] - x[None, :]) / width) ** 2)
    return (w @ y) / w.sum(axis=1)


def make_disk(gp: GalaxyParams = GalaxyParams(), ic: DiskICParams = DiskICParams(),
              num: NumericalParams = NumericalParams(), verbose=True) -> DiskParticles:
    """Öz-yerçekimli disk parçacıklarını üretir. Kütle merkezi ve net momentum sıfırlanır.

    verbose=True iken ELN çubuk kararlılık tanısını yazdırır.
    """
    if verbose:
        e = eln_parameter(gp, ic)
        verdict = "çubuk oluşması beklenir" if e < ELN_BAR_THRESHOLD else "çubuğa karşı kararlı"
        print(f"[IC] M_d = {ic.M_d:.4e} M_sun, Q_target = {ic.Q_target}, N = {num.N}; "
              f"ELN eps = V_max/sqrt(G M_d/R_d) = {e:.3f} (< {ELN_BAR_THRESHOLD}: {verdict})")
    rng = np.random.default_rng(num.seed)
    N = num.N

    R = sample_radii(rng, N, ic.R_d, ic.R_cut)
    z = sample_heights(rng, N, ic.z_0)
    phi = rng.uniform(0.0, 2.0 * np.pi, N)

    c, s = np.cos(phi), np.sin(phi)
    pos = np.column_stack((R * c, R * s, z))
    mass = np.full(N, ic.M_d / N)

    R_grid = np.geomspace(1e-3 * ic.R_d, ic.R_cut, num.n_profile_grid)
    prof = disk_profiles(gp, ic, R_grid)
    at = lambda key: np.interp(R, R_grid, prof[key])

    # Gerçek disk kalın (sech^2) ve yumuşatılmış: düzlemdeki v_c, ince Freeman diskinden
    # düşüktür. <v_phi>, parçacıkların kendi alanından ölçülen v_c ile kurulur.
    R_ring = np.linspace(ic.R_cut / num.ic_ring_n, ic.R_cut, num.ic_ring_n)
    vc2_an = circular_velocity(ic_potential(gp, ic), R_ring) ** 2
    dvc2 = _gaussian_smooth(R_ring, particle_vc2(gp, num, pos, mass, R_ring) - vc2_an,
                            num.ic_vc_smooth)
    vc2 = at("v_c") ** 2 + np.interp(R, R_ring, dvc2)
    vphi_mean = np.sqrt(np.maximum(vc2 + at("sigma_R") ** 2 * at("drift_term"), 0.0))

    v_R = rng.normal(0.0, 1.0, N) * at("sigma_R")
    v_phi = vphi_mean + rng.normal(0.0, 1.0, N) * at("sigma_phi")
    v_z = rng.normal(0.0, 1.0, N) * at("sigma_z")
    vel = np.column_stack((v_R * c - v_phi * s, v_R * s + v_phi * c, v_z))

    pos -= np.average(pos, axis=0, weights=mass)
    vel -= np.average(vel, axis=0, weights=mass)
    return DiskParticles(pos, vel, mass, R_ring, dvc2)


def measure_profiles(p: DiskParticles, edges):
    """Parçacıklardan halka bazında ölçüm. edges : kpc (halka sınırları).

    dönüş : dict — R_mid [kpc], Sigma [M_sun/kpc^2], sigma_R, vphi_mean [kpc/Myr], count
    """
    x, y = p.pos[:, 0], p.pos[:, 1]
    R = np.hypot(x, y)
    v_R = (x * p.vel[:, 0] + y * p.vel[:, 1]) / R
    v_phi = (x * p.vel[:, 1] - y * p.vel[:, 0]) / R

    idx = np.digitize(R, edges) - 1
    nb = len(edges) - 1
    ok = (idx >= 0) & (idx < nb)
    count = np.bincount(idx[ok], minlength=nb)
    m = np.bincount(idx[ok], weights=p.mass[ok], minlength=nb)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)

    def mean(v):
        return np.bincount(idx[ok], weights=v[ok], minlength=nb) / np.maximum(count, 1)

    mR = mean(v_R)
    sigR = np.sqrt(np.maximum(mean(v_R ** 2) - mR ** 2, 0.0))
    return dict(R_mid=0.5 * (edges[1:] + edges[:-1]), Sigma=m / area, sigma_R=sigR,
                vphi_mean=mean(v_phi), count=count)


def asymmetric_drift_check(gp: GalaxyParams, ic: DiskICParams, disk: DiskParticles, R0, half_width):
    """v_c - <v_phi> farkı: parçacıklardan ölçülen ve Jeans denkleminin öngördüğü.

    R0 - half_width < R < R0 + half_width halkasındaki parçacıklar üzerinden:
      v_c       : parçacıklardan ölçülen (analitik + disk.vc2_correction)
      ölçülen   : < v_c(R_i) - v_phi,i >
      öngörülen : < v_c - sqrt(v_c^2 + sigma_R^2 D) >,  D = dln(Sigma sigma_R^2)/dlnR + 1 - kappa^2/4Omega^2
      1. mertebe: < sigma_R^2 (-D) / (2 v_c) >
    R0, half_width : kpc  ->  dict, hızlar km/s
    """
    x, y = disk.pos[:, 0], disk.pos[:, 1]
    R = np.hypot(x, y)
    sel = np.abs(R - R0) < half_width
    R, vphi = R[sel], ((x * disk.vel[:, 1] - y * disk.vel[:, 0]) / np.hypot(x, y))[sel]
    prof = disk_profiles(gp, ic, R)
    vc = np.sqrt(prof["v_c"] ** 2 + disk.vc2_correction(R))
    s2D = prof["sigma_R"] ** 2 * prof["drift_term"]
    return dict(
        n=int(sel.sum()),
        v_c=kpcmyr_to_kms(vc.mean()),
        vphi_mean=kpcmyr_to_kms(vphi.mean()),
        measured=kpcmyr_to_kms(np.mean(vc - vphi)),
        measured_err=kpcmyr_to_kms(vphi.std() / np.sqrt(sel.sum())),
        predicted=kpcmyr_to_kms(np.mean(vc - np.sqrt(vc ** 2 + s2D))),
        predicted_first_order=kpcmyr_to_kms(np.mean(-s2D / (2.0 * vc))),
    )


def measured_Q(gp: GalaxyParams, ic: DiskICParams, meas):
    """Ölçülen sigma_R ve Sigma ile, analitik kappa kullanarak Q(R)."""
    ka = kappa(ic_potential(gp, ic), meas["R_mid"])
    return meas["sigma_R"] * ka / (TOOMRE_STELLAR * G * meas["Sigma"])
