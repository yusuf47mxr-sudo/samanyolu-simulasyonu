"""Sabit analitik potansiyelde test parçacık yörüngeleri — leapfrog (KDK).

Leapfrog simplektiktir: enerji hatası sınırlı kalır, sürüklenmez.
Eksenel simetrik potansiyelde kick yalnızca radyal (R yönü) ve z yönünde
kuvvet uygular, drift de L_z'yi değiştirmez; bu yüzden L_z yuvarlama
hatası düzeyinde korunur.

Birimler: konum : kpc, hız : kpc/Myr, zaman : Myr, enerji : kpc^2/Myr^2 (birim kütle başına).
Diziler (N, 3) biçimindedir: sütunlar x, y, z.
"""

from dataclasses import dataclass

import numpy as np

from .rotation import circular_velocity
from .units import kms_to_kpcmyr


def acceleration(pot, pos):
    """a = -grad Phi.  pos : (N,3) kpc  ->  (N,3) kpc/Myr^2."""
    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    R = np.hypot(x, y)
    fR = -pot.dphi_dR(R, z) / R
    return np.column_stack((fR * x, fR * y, -pot.dphi_dz(R, z)))


def energy(pot, pos, vel):
    """Birim kütle başına toplam enerji. pos : kpc, vel : kpc/Myr -> kpc^2/Myr^2."""
    R = np.hypot(pos[..., 0], pos[..., 1])
    return 0.5 * np.sum(vel * vel, axis=-1) + pot.phi(R, pos[..., 2])


def angular_momentum_z(pos, vel):
    """L_z = x v_y - y v_x.  pos : kpc, vel : kpc/Myr -> kpc^2/Myr."""
    return pos[..., 0] * vel[..., 1] - pos[..., 1] * vel[..., 0]


@dataclass
class Orbit:
    """Leapfrog çıktısı. t : (M,) Myr; pos, vel : (M, N, 3) kpc, kpc/Myr."""
    t: np.ndarray
    pos: np.ndarray
    vel: np.ndarray


def leapfrog(pot, pos0, vel0, dt, n_steps, save_every=1):
    """Kick–drift–kick leapfrog.

    pot        : potansiyel nesnesi
    pos0, vel0 : (N,3) ya da (3,) — kpc, kpc/Myr
    dt         : Myr
    n_steps    : adım sayısı
    save_every : her kaç adımda bir kayıt alınacağı (ilk ve son durum her zaman kaydedilir)
    dönüş      : Orbit
    """
    pos = np.atleast_2d(np.asarray(pos0, dtype=float)).copy()
    vel = np.atleast_2d(np.asarray(vel0, dtype=float)).copy()

    n_save = n_steps // save_every + 1 + (n_steps % save_every != 0)
    t_out = np.empty(n_save)
    pos_out = np.empty((n_save,) + pos.shape)
    vel_out = np.empty((n_save,) + vel.shape)
    t_out[0], pos_out[0], vel_out[0] = 0.0, pos, vel
    k = 1

    acc = acceleration(pot, pos)
    for i in range(1, n_steps + 1):
        vel += 0.5 * dt * acc
        pos += dt * vel
        acc = acceleration(pot, pos)
        vel += 0.5 * dt * acc
        if i % save_every == 0 or i == n_steps:
            t_out[k], pos_out[k], vel_out[k] = i * dt, pos, vel
            k += 1
    return Orbit(t_out[:k], pos_out[:k], vel_out[:k])


def initial_conditions(pot, R0, v_phi_frac, v_R0_kms=0.0, v_z0_kms=0.0):
    """Disk düzleminde (x=R0, y=0, z=0) başlayan tek parçacık.

    R0 : kpc; v_phi_frac : v_phi / v_c(R0); v_R0_kms, v_z0_kms : km/s
    dönüş : pos (1,3) kpc, vel (1,3) kpc/Myr
    """
    v_phi = v_phi_frac * float(circular_velocity(pot, R0))
    pos = np.array([[R0, 0.0, 0.0]])
    vel = np.array([[kms_to_kpcmyr(v_R0_kms), v_phi, kms_to_kpcmyr(v_z0_kms)]])
    return pos, vel
