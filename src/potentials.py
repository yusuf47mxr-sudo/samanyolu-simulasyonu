"""Analitik galaksi potansiyelleri: Miyamoto–Nagai, Hernquist, NFW, Kepler.

Tüm fonksiyonlar silindirik (R, z) koordinatında çalışır ve numpy dizilerini
kabul eder. Birimler: R, z, r : kpc; Phi : kpc^2/Myr^2; kütle : M_sun.

Her bileşen şunları sağlar:
    phi(R, z)               potansiyel
    dphi_dR(R, z)           dPhi/dR   (analitik)
    dphi_dz(R, z)           dPhi/dz   (analitik)
    d2phi_dR2_midplane(R)   d2Phi/dR2 at z=0 (analitik; kappa için)
    enclosed_mass(r)        r yarıçaplı küre içindeki kütle
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import ive, kve

from .config import GalaxyParams
from .units import G


class MiyamotoNagai:
    """Miyamoto–Nagai disk. M : M_sun, a : kpc, b : kpc."""

    name = "disk"

    def __init__(self, M, a, b):
        self.M, self.a, self.b = M, a, b
        self.GM = G * M

    def _zeta(self, z):
        return np.sqrt(z * z + self.b * self.b)

    def _D2(self, R, z):
        return R * R + (self.a + self._zeta(z)) ** 2

    def phi(self, R, z):
        """R, z : kpc -> Phi : kpc^2/Myr^2."""
        return -self.GM / np.sqrt(self._D2(R, z))

    def dphi_dR(self, R, z):
        """R, z : kpc -> dPhi/dR : kpc/Myr^2."""
        return self.GM * R / self._D2(R, z) ** 1.5

    def dphi_dz(self, R, z):
        """R, z : kpc -> dPhi/dz : kpc/Myr^2."""
        zeta = self._zeta(z)
        return self.GM * z * (self.a + zeta) / (zeta * self._D2(R, z) ** 1.5)

    def density(self, R, z):
        """R, z : kpc -> rho : M_sun/kpc^3 (Poisson denkleminden analitik)."""
        a, b = self.a, self.b
        zeta = self._zeta(z)
        num = a * R * R + (a + 3.0 * zeta) * (a + zeta) ** 2
        den = self._D2(R, z) ** 2.5 * zeta ** 3
        return b * b * self.M / (4.0 * np.pi) * num / den

    def d2phi_dR2_midplane(self, R):
        """R : kpc -> d2Phi/dR2 (z=0) : 1/Myr^2."""
        s2 = (self.a + self.b) ** 2
        return self.GM * (s2 - 2.0 * R * R) / (R * R + s2) ** 2.5

    def enclosed_mass(self, r):
        """r : kpc -> küre içindeki kütle : M_sun.

        Disk küresel değildir; Gauss yasasıyla küre yüzeyinden geçen akı
        integrali alınır:  M = r^2/G * int_0^{pi/2} (dPhi/dr) sin(th) dth.
        Bu, rho'nun küre hacmi üzerindeki integraline eşittir; ince diskte
        hacim integralinden çok daha kararlıdır (testte karşılaştırılır).
        """
        def one(rr):
            def integrand(th):
                R, z = rr * np.sin(th), rr * np.cos(th)
                dphi_dr = self.dphi_dR(R, z) * np.sin(th) + self.dphi_dz(R, z) * np.cos(th)
                return dphi_dr * np.sin(th)
            val, _ = quad(integrand, 0.0, 0.5 * np.pi, limit=200)
            return rr * rr * val / G

        return np.vectorize(one, otypes=[float])(r)


class FreemanDisk:
    """İnce (sonsuz ince) üstel disk, Freeman (1970). Yalnızca z = 0 düzlemi.

    Sigma(R) = Sigma_0 exp(-R/R_d). Sigma_0, R_cut içindeki kütle M olacak
    şekilde seçilir (kesilmiş N-cisim diskiyle aynı normalizasyon); potansiyel
    ise kesilmemiş diskinkidir — R_cut = 15 kpc'de dışarıda kalan kütle %1.7.

    M : M_sun, R_d, R_cut : kpc. Başlangıç koşulları için kullanılır;
    dphi_dz tanımlı değildir, yörünge entegrasyonunda kullanılmaz.
    """

    name = "disk"

    def __init__(self, M, R_d, R_cut):
        self.M, self.R_d, self.R_cut = M, R_d, R_cut
        x = R_cut / R_d
        self.Sigma0 = M / (2.0 * np.pi * R_d ** 2 * (1.0 - (1.0 + x) * np.exp(-x)))
        self.GS = G * self.Sigma0

    @staticmethod
    def _bessel(y):
        # I_n(y) K_m(y) çarpımları; ölçekli fonksiyonlarla taşma olmaz
        i0, i1 = ive(0, y), ive(1, y)
        k0, k1 = kve(0, y), kve(1, y)
        return i0, i1, k0, k1

    def phi(self, R, z=0.0):
        """R : kpc (z=0) -> Phi : kpc^2/Myr^2."""
        y = np.asarray(R, dtype=float) / (2.0 * self.R_d)
        i0, i1, k0, k1 = self._bessel(y)
        return -np.pi * self.GS * R * (i0 * k1 - i1 * k0)

    def _vc2(self, y):
        i0, i1, k0, k1 = self._bessel(y)
        return 4.0 * np.pi * self.GS * self.R_d * y * y * (i0 * k0 - i1 * k1)

    def dphi_dR(self, R, z=0.0):
        """R : kpc (z=0) -> dPhi/dR : kpc/Myr^2."""
        R = np.asarray(R, dtype=float)
        return self._vc2(R / (2.0 * self.R_d)) / R

    def dphi_dz(self, R, z):
        raise NotImplementedError("FreemanDisk yalnızca z = 0 düzleminde tanımlı")

    def d2phi_dR2_midplane(self, R):
        """R : kpc -> d2Phi/dR2 (z=0) : 1/Myr^2."""
        R = np.asarray(R, dtype=float)
        y = R / (2.0 * self.R_d)
        i0, i1, k0, k1 = self._bessel(y)
        g = i0 * k0 - i1 * k1
        dg = 2.0 * (i1 * k0 - i0 * k1 + i1 * k1 / y)
        dvc2_dR = 2.0 * np.pi * self.GS * (2.0 * y * g + y * y * dg)
        return dvc2_dR / R - self._vc2(y) / (R * R)

    def surface_density(self, R):
        """R : kpc -> Sigma : M_sun/kpc^2 (R > R_cut için 0)."""
        R = np.asarray(R, dtype=float)
        return np.where(R <= self.R_cut, self.Sigma0 * np.exp(-R / self.R_d), 0.0)


class _Spherical:
    """Küresel potansiyeller için ortak (R, z) arayüzü."""

    def _dphi_dr(self, r):
        raise NotImplementedError

    def dphi_dR(self, R, z):
        """R, z : kpc -> dPhi/dR : kpc/Myr^2."""
        r = np.hypot(R, z)
        return self._dphi_dr(r) * R / r

    def dphi_dz(self, R, z):
        """R, z : kpc -> dPhi/dz : kpc/Myr^2."""
        r = np.hypot(R, z)
        return self._dphi_dr(r) * z / r


class Hernquist(_Spherical):
    """Hernquist şişkinlik. M : M_sun, c : kpc."""

    name = "bulge"

    def __init__(self, M, c):
        self.M, self.c = M, c
        self.GM = G * M

    def phi(self, R, z):
        """R, z : kpc -> Phi : kpc^2/Myr^2."""
        return -self.GM / (np.hypot(R, z) + self.c)

    def _dphi_dr(self, r):
        return self.GM / (r + self.c) ** 2

    def d2phi_dR2_midplane(self, R):
        """R : kpc -> d2Phi/dR2 (z=0) : 1/Myr^2."""
        return -2.0 * self.GM / (R + self.c) ** 3

    def enclosed_mass(self, r):
        """r : kpc -> M(<r) : M_sun."""
        return self.M * r * r / (r + self.c) ** 2


class NFW(_Spherical):
    """NFW karanlık madde halesi. M_vir : M_sun, r_s : kpc, conc = r_vir/r_s."""

    name = "halo"

    def __init__(self, M_vir, r_s, conc):
        self.M_vir, self.r_s, self.conc = M_vir, r_s, conc
        # 4*pi*rho_s*r_s^3 = M_vir / f(conc)   [M_sun]
        self.M_norm = M_vir / self._f(conc)
        self.GM = G * self.M_norm

    @staticmethod
    def _f(x):
        """ln(1+x) - x/(1+x), boyutsuz."""
        return np.log1p(x) - x / (1.0 + x)

    @property
    def rho_s(self):
        """Karakteristik yoğunluk : M_sun/kpc^3."""
        return self.M_norm / (4.0 * np.pi * self.r_s ** 3)

    def phi(self, R, z):
        """R, z : kpc -> Phi : kpc^2/Myr^2."""
        r = np.hypot(R, z)
        return -self.GM * np.log1p(r / self.r_s) / r

    def _dphi_dr(self, r):
        return self.GM * self._f(r / self.r_s) / (r * r)

    def d2phi_dR2_midplane(self, R):
        """R : kpc -> d2Phi/dR2 (z=0) : 1/Myr^2."""
        x = R / self.r_s
        return self.GM * (x / ((1.0 + x) ** 2 * self.r_s * R * R) - 2.0 * self._f(x) / R ** 3)

    def enclosed_mass(self, r):
        """r : kpc -> M(<r) : M_sun."""
        return self.M_norm * self._f(r / self.r_s)


class Kepler(_Spherical):
    """Nokta kütle (merkezi kara delik). M : M_sun, eps : kpc.

    eps > 0 ise Plummer yumuşatması: Phi = -GM / sqrt(r^2 + eps^2).
    Analitik model eps = 0 kullanır; N-cisimde parçacıklarla aynı eps verilir.
    """

    name = "bh"

    def __init__(self, M, eps=0.0):
        self.M, self.eps = M, eps
        self.GM = G * M

    def phi(self, R, z):
        """R, z : kpc -> Phi : kpc^2/Myr^2."""
        return -self.GM / np.sqrt(R * R + z * z + self.eps ** 2)

    def _dphi_dr(self, r):
        return self.GM * r / (r * r + self.eps ** 2) ** 1.5

    def d2phi_dR2_midplane(self, R):
        """R : kpc -> d2Phi/dR2 (z=0) : 1/Myr^2."""
        e2 = self.eps ** 2
        return self.GM * (e2 - 2.0 * R * R) / (R * R + e2) ** 2.5

    def enclosed_mass(self, r):
        """r : kpc -> M(<r) : M_sun (nokta kütle; yumuşatma kütleyi dağıtmaz)."""
        return self.M * np.ones_like(np.asarray(r, dtype=float))


class CompositePotential:
    """Bileşenlerin toplamı. Potansiyel ve kuvvetler toplanır."""

    def __init__(self, components):
        self.components = list(components)

    def __getitem__(self, name):
        for c in self.components:
            if c.name == name:
                return c
        raise KeyError(name)

    def names(self):
        return [c.name for c in self.components]

    def phi(self, R, z):
        """R, z : kpc -> Phi : kpc^2/Myr^2."""
        return sum(c.phi(R, z) for c in self.components)

    def dphi_dR(self, R, z):
        """R, z : kpc -> dPhi/dR : kpc/Myr^2."""
        return sum(c.dphi_dR(R, z) for c in self.components)

    def dphi_dz(self, R, z):
        """R, z : kpc -> dPhi/dz : kpc/Myr^2."""
        return sum(c.dphi_dz(R, z) for c in self.components)

    def d2phi_dR2_midplane(self, R):
        """R : kpc -> d2Phi/dR2 (z=0) : 1/Myr^2."""
        return sum(c.d2phi_dR2_midplane(R) for c in self.components)

    def enclosed_mass(self, r):
        """r : kpc -> küre içindeki toplam kütle : M_sun."""
        return sum(c.enclosed_mass(r) for c in self.components)


def build_potential(p: GalaxyParams = GalaxyParams(), eps_bh=0.0) -> CompositePotential:
    """Fiziksel parametrelerden toplam potansiyeli kurar.

    Hale kapalıysa (p.halo.enabled=False) NFW bileşeni hiç eklenmez.
    eps_bh : kpc, kara delik yumuşatması (analitik modelde 0).
    """
    comps = [
        Kepler(p.bh.M_bh, eps_bh),
        Hernquist(p.bulge.M_b, p.bulge.c),
        MiyamotoNagai(p.disk.M_d, p.disk.a, p.disk.b),
    ]
    if p.halo.enabled:
        comps.append(NFW(p.halo.M_vir, p.halo.r_s, p.halo.conc))
    return CompositePotential(comps)


def build_external_potential(p: GalaxyParams = GalaxyParams(), eps_bh=0.0) -> CompositePotential:
    """N-cisimde sabit kalan kısım: kara delik + şişkinlik + (açıksa) hale.

    Disk dahil değildir; disk parçacıklardan öz-yerçekimiyle hesaplanır.
    eps_bh : kpc, kara delik yumuşatması.
    """
    full = build_potential(p, eps_bh)
    return CompositePotential(c for c in full.components if c.name != "disk")
