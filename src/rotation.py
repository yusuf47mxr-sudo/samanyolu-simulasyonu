"""Dönme eğrisi, açısal frekans Omega ve epicycle frekansı kappa.

Tüm türevler analitiktir (potentials.py). Sayısal türev yalnızca
test karşılaştırması için `circular_velocity_numerical` içinde durur.

Birimler: R : kpc, v : kpc/Myr, Omega, kappa : 1/Myr, kütle : M_sun.
"""

from dataclasses import replace

import numpy as np
from scipy.optimize import brentq

from .config import HaloParams, param_range
from .potentials import build_potential
from .units import G, kms_to_kpcmyr


def circular_velocity(pot, R):
    """v_c(R) = sqrt(R * dPhi/dR |_{z=0}).

    pot : potansiyel nesnesi (bileşen ya da CompositePotential)
    R   : kpc  ->  v_c : kpc/Myr
    """
    R = np.asarray(R, dtype=float)
    return np.sqrt(np.maximum(R * pot.dphi_dR(R, 0.0), 0.0))


def circular_velocity_components(pot, R):
    """Bileşen bazında v_c. R : kpc -> {ad: v_c [kpc/Myr]}."""
    return {c.name: circular_velocity(c, R) for c in pot.components}


def circular_velocity_numerical(pot, R, h=1e-4):
    """Merkezi fark ile v_c — yalnızca test amaçlı. R, h : kpc -> kpc/Myr."""
    R = np.asarray(R, dtype=float)
    dphi = (pot.phi(R + h, 0.0) - pot.phi(R - h, 0.0)) / (2.0 * h)
    return np.sqrt(R * dphi)


def omega(pot, R):
    """Omega(R) = v_c / R. R : kpc -> 1/Myr."""
    return circular_velocity(pot, R) / np.asarray(R, dtype=float)


def kappa(pot, R):
    """Epicycle frekansı. R : kpc -> kappa : 1/Myr.

    kappa^2 = R dOmega^2/dR + 4 Omega^2 = Phi'' + 3 Phi'/R   (z=0)
    """
    R = np.asarray(R, dtype=float)
    k2 = pot.d2phi_dR2_midplane(R) + 3.0 * pot.dphi_dR(R, 0.0) / R
    return np.sqrt(np.maximum(k2, 0.0))


def orbital_period(pot, R):
    """Dairesel yörünge periyodu 2 pi R / v_c. R : kpc -> Myr."""
    R = np.asarray(R, dtype=float)
    return 2.0 * np.pi * R / circular_velocity(pot, R)


def escape_velocity(pot, R, z=0.0):
    """v_esc = sqrt(-2 Phi). R, z : kpc -> kpc/Myr."""
    return np.sqrt(np.maximum(-2.0 * pot.phi(R, z), 0.0))


def dynamical_mass(pot, R):
    """Küresel yaklaşımla dinamik kütle R v_c^2 / G. R : kpc -> M_sun.

    Gözlemsel 'M(<R)' değerleri genelde bu tahminle verilir.
    """
    R = np.asarray(R, dtype=float)
    return R * circular_velocity(pot, R) ** 2 / G


def solve_mvir_for_vc(params, v_target_kms, R, M_vir_range=None):
    """Ters mod: v_c(R) = v_target olacak M_vir'i brentq ile çözer.

    params        : GalaxyParams (hale açık olmalı)
    v_target_kms  : km/s  (arayüz girdisi)
    R             : kpc
    M_vir_range   : (alt, üst) M_sun; verilmezse config'deki aralık
    dönüş         : M_vir [M_sun]; aralıkta çözüm yoksa None
    """
    lo, hi = M_vir_range or param_range(HaloParams, "M_vir")
    v_t = kms_to_kpcmyr(v_target_kms)
    halo_on = replace(params.halo, enabled=True)

    def resid(M_vir):
        p = replace(params, halo=replace(halo_on, M_vir=M_vir))
        return float(circular_velocity(build_potential(p), R)) - v_t

    f_lo, f_hi = resid(lo), resid(hi)
    if f_lo * f_hi > 0:
        return None
    return brentq(resid, lo, hi, xtol=1e-6 * hi, rtol=1e-12)


def exponential_surface_density(R, M_disk, R_d, R_cut):
    """Kesilmiş üstel disk Sigma(R) = Sigma_0 exp(-R/R_d), R < R_cut.

    R, R_d, R_cut : kpc;  M_disk : M_sun  ->  Sigma : M_sun/kpc^2
    Sigma_0, R_cut içindeki toplam kütle M_disk olacak şekilde seçilir.
    """
    x = R_cut / R_d
    sigma0 = M_disk / (2.0 * np.pi * R_d ** 2 * (1.0 - (1.0 + x) * np.exp(-x)))
    R = np.asarray(R, dtype=float)
    return np.where(R < R_cut, sigma0 * np.exp(-R / R_d), 0.0)


def toomre_q(pot, R, sigma_R, Sigma):
    """Q = sigma_R kappa / (3.36 G Sigma).

    R : kpc, sigma_R : kpc/Myr, Sigma : M_sun/kpc^2  ->  Q : boyutsuz
    """
    return sigma_R * kappa(pot, R) / (3.36 * G * Sigma)


def sigma_R_profile(pot, R, Sigma_fn, Q_target, R_ref, R_d):
    """sigma_R(R) ∝ exp(-R / 2R_d), Q(R_ref) = Q_target olacak şekilde ölçeklenir.

    R, R_ref, R_d : kpc; Sigma_fn : R -> M_sun/kpc^2  ->  sigma_R : kpc/Myr
    """
    s_ref = Q_target * 3.36 * G * Sigma_fn(R_ref) / kappa(pot, R_ref)
    return s_ref * np.exp(-(np.asarray(R, dtype=float) - R_ref) / (2.0 * R_d))
