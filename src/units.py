"""Birim sistemi ve dönüşümler.

İç birimler: uzunluk kpc, zaman Myr, kütle M_sun.
Hızlar içeride kpc/Myr'dir; km/s'ye yalnızca ekrana yazarken çevrilir.
"""

G = 4.498502e-12        # kpc^3 / (M_sun * Myr^2)
KMS_KPMY = 1.022712e-3  # 1 km/s -> kpc/Myr


def kms_to_kpcmyr(v):
    """km/s -> kpc/Myr.

    v : km/s  ->  dönüş : kpc/Myr
    """
    return v * KMS_KPMY


def kpcmyr_to_kms(v):
    """kpc/Myr -> km/s.

    v : kpc/Myr  ->  dönüş : km/s
    """
    return v / KMS_KPMY
