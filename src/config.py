"""Tüm varsayılan parametreler — tek yer.

Fiziksel ve sayısal parametreler ayrı dataclass'larda durur (PROJE_TANIMI.md §8).
Her alanın `metadata` içinde birimi ve izinli aralığı vardır; arayüz
kaydıraklarını bu aralıklardan kurar.

Birimler: uzunluk kpc, zaman Myr, kütle M_sun.
"""

from dataclasses import dataclass, field, fields


def _p(default, lo, hi, unit):
    """Aralık ve birim bilgisi taşıyan dataclass alanı."""
    return field(default=default, metadata={"range": (lo, hi), "unit": unit})


def param_range(cls, name):
    """Bir parametrenin (alt, üst) aralığı."""
    return next(f.metadata["range"] for f in fields(cls) if f.name == name)


def param_unit(cls, name):
    """Bir parametrenin birimi."""
    return next(f.metadata["unit"] for f in fields(cls) if f.name == name)


# --------------------------------------------------------------------------
# Fiziksel parametreler
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class DiskParams:
    """Miyamoto–Nagai disk. M_d : M_sun, a : kpc, b : kpc."""
    M_d: float = _p(6.8e10, 3e10, 12e10, "M_sun")
    a: float = _p(3.0, 1.5, 6.0, "kpc")
    b: float = _p(0.28, 0.1, 1.0, "kpc")


@dataclass(frozen=True)
class BulgeParams:
    """Hernquist şişkinlik. M_b : M_sun, c : kpc."""
    M_b: float = _p(5.0e9, 0.0, 2e10, "M_sun")
    c: float = _p(0.5, 0.1, 2.0, "kpc")


@dataclass(frozen=True)
class HaloParams:
    """NFW karanlık madde halesi. M_vir : M_sun, r_s : kpc, conc : boyutsuz.

    enabled=False iken hale potansiyele hiç katılmaz.
    """
    M_vir: float = _p(1.0e12, 0.0, 3e12, "M_sun")
    r_s: float = _p(16.0, 8.0, 30.0, "kpc")
    conc: float = _p(15.0, 5.0, 25.0, "")
    enabled: bool = True


@dataclass(frozen=True)
class BlackHoleParams:
    """Merkezi kara delik (Kepler). M_bh : M_sun."""
    M_bh: float = _p(4.3e6, 0.0, 1e10, "M_sun")


@dataclass(frozen=True)
class GalaxyParams:
    """Toplam potansiyeli tanımlayan fiziksel parametreler."""
    disk: DiskParams = field(default_factory=DiskParams)
    bulge: BulgeParams = field(default_factory=BulgeParams)
    halo: HaloParams = field(default_factory=HaloParams)
    bh: BlackHoleParams = field(default_factory=BlackHoleParams)


@dataclass(frozen=True)
class DiskICParams:
    """N-cisim disk başlangıç koşulları (fiziksel). R_d, R_cut, z_0 : kpc; t_end : Myr."""
    # Üstel disk kütlesi. Varsayılan sabit seçilmedi: ics.solve_disk_mass ile
    # üstel disk + şişkinlik + hale + kara delik toplamı v_c(R_SUN) = VC_SUN_TARGET
    # verecek şekilde brentq ile çözüldü (tests/test_ics.py bunu doğrular).
    # MN'nin M_d = 6.8e10'u bir uyum fonksiyonu parametresiydi; burada geçerli değil.
    M_d: float = _p(4.82486e10, 3e10, 12e10, "M_sun")
    Q_target: float = _p(1.3, 0.8, 2.5, "")
    R_d: float = _p(2.6, 2.6, 2.6, "kpc")
    R_cut: float = _p(15.0, 15.0, 15.0, "kpc")
    z_0: float = _p(0.3, 0.3, 0.3, "kpc")
    t_end: float = _p(2000.0, 500.0, 3000.0, "Myr")


# --------------------------------------------------------------------------
# Sayısal parametreler — sonucu değiştirmemelidir
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class NumericalParams:
    """N : parçacık sayısı, eps : kpc, dt : Myr, theta : Barnes–Hut açılma açısı."""
    N: int = _p(20000, 5000, 50000, "")
    eps: float = _p(0.08, 0.05, 0.15, "kpc")
    dt: float = _p(0.2, 0.05, 0.5, "Myr")
    theta: float = _p(0.6, 0.4, 0.8, "")
    seed: int = 42                      # rastgele sayı üreteci tohumu
    snap_interval: float = 10.0         # Myr, snapshot yazma aralığı
    energy_interval: float = 5.0        # Myr, enerji tanısı aralığı (doğrudan toplam, O(N^2))
    leaf_size: int = 8                  # Barnes–Hut yaprak kovası kapasitesi
    tree_quadrupole: bool = True        # ağaçta kuadrupol terimi (theta=0.6'da RMS kuvvet hatası ~%0.1)
    n_profile_grid: int = 2000          # IC radyal profil ızgarası nokta sayısı
    ic_ring_n: int = 120                # IC: parçacık v_c'sinin ölçüldüğü halka sayısı
    ic_ring_nphi: int = 64              # IC: halka başına azimut noktası
    ic_vc_smooth: float = 0.5           # kpc, IC v_c^2 düzeltmesinin Gauss düzgünleştirme genişliği


# --------------------------------------------------------------------------
# Gözlemsel referanslar
# --------------------------------------------------------------------------

R_SUN = 8.2  # kpc, Güneş'in galaktosentrik uzaklığı
VC_SUN_TARGET = 228.7  # km/s, N-cisim disk kütlesinin normalize edildiği v_c(R_SUN)
Q_REF_RADIUS = R_SUN  # kpc, Q(R) profilinde Q = Q_target olan referans yarıçap


@dataclass(frozen=True)
class OrbitParams:
    """Tek test parçacık yörüngesi (fiziksel).

    R0 : kpc, v_phi_frac : v_phi / v_c(R0), v_R0 ve v_z0 : km/s, t_orbit : Myr.
    """
    R0: float = _p(8.2, 1.0, 30.0, "kpc")
    v_phi_frac: float = _p(0.8, 0.0, 1.5, "")
    v_R0: float = _p(30.0, -200.0, 200.0, "km/s")
    v_z0: float = _p(40.0, -200.0, 200.0, "km/s")
    t_orbit: float = _p(1000.0, 100.0, 5000.0, "Myr")


# Çubuk gücü A2/A0'ın ölçüldüğü bölge (kpc); N-cisim koşusu sırasında hesaplanır
A2_R_MAX = 5.0

# Ters mod varsayılan hedefi
V_TARGET_DEFAULT = 229.0  # km/s, v_c(R_SUN) hedefi

# Gözlem: Eilers+2019 (ApJ 871, 120) doğrusal uyumu
#   v_c(R) = 229.0 - 1.7 (R - 8.122) km/s,  5 <= R <= 25 kpc
# Noktalar bu uyumdan örneklenir; tek tek veri noktası değildir.
OBS_EILERS2019 = dict(v0=229.0, slope=-1.7, R0=8.122, R_min=5.0, R_max=25.0, n=11, err=5.0)

# Arayüzde çizilen yarıçap aralıkları (kpc)
PLOT_R_LIN = (0.1, 30.0)
PLOT_R_LOG = (1e-3, 300.0)
