"""Samanyolu — etkileşimli dönme eğrisi ve yörünge arayüzü (Aşama 1).

Çalıştırma:  .venv/bin/streamlit run app.py

İç hesaplar kpc / Myr / M_sun birimlerindedir; km/s'ye yalnızca burada,
ekrana yazarken çevrilir.
"""

import base64
import os
import json
from dataclasses import fields
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

from src import analysis as an
from src.config import (
    A2_R_MAX,
    OBS_EILERS2019,
    PLOT_R_LIN,
    PLOT_R_LOG,
    Q_REF_RADIUS,
    R_SUN,
    V_TARGET_DEFAULT,
    BlackHoleParams,
    BulgeParams,
    DiskICParams,
    DiskParams,
    GalaxyParams,
    HaloParams,
    NumericalParams,
    OrbitParams,
    param_range,
)
from src.ics import asymmetric_drift_check, make_disk, sample_heights, sample_radii
from src.orbits import angular_momentum_z, energy, initial_conditions, leapfrog
from src.potentials import build_potential
from src.rotation import (
    circular_velocity,
    circular_velocity_components,
    dynamical_mass,
    escape_velocity,
    exponential_surface_density,
    orbital_period,
    sigma_R_profile,
    solve_mvir_for_vc,
    toomre_q,
)
from src.units import kpcmyr_to_kms

# Arayüz görünümü (fizik değil): bileşen renkleri ve etiketleri
COLORS = {"total": "#111111", "bh": "#7b3294", "bulge": "#d95f02", "disk": "#1b9e77", "halo": "#386cb0"}
LABELS = {"total": "Toplam", "bh": "Kara delik", "bulge": "Şişkinlik", "disk": "Disk", "halo": "Karanlık madde halesi"}
MAX_ORBIT_SAMPLES = 4000  # yörünge grafiğinde saklanan en fazla nokta
DRIFT_HALF_WIDTH = 1.0    # kpc, sürüklenme ölçüm halkasının yarı genişliği

st.set_page_config(page_title="Samanyolu Dinamiği", layout="wide")


def default(cls, name):
    return next(f.default for f in fields(cls) if f.name == name)


def mass_slider(container, label, cls, name, scale, unit_label, step, **kw):
    """M_sun cinsinden parametre için ölçekli kaydırak; dönüş M_sun."""
    lo, hi = param_range(cls, name)
    v = container.slider(
        f"{label} [{unit_label}]", lo / scale, hi / scale, default(cls, name) / scale, step, **kw
    )
    return v * scale


def plain_slider(container, label, cls, name, step, unit=None, **kw):
    lo, hi = param_range(cls, name)
    lab = f"{label} [{unit}]" if unit else label
    return container.slider(lab, float(lo), float(hi), float(default(cls, name)), step, **kw)


# ==========================================================================
# Kenar çubuğu — kontroller
# ==========================================================================

st.sidebar.title("Parametreler")
tab_main, tab_adv, tab_num = st.sidebar.tabs(["Ana", "Gelişmiş", "Sayısal"])

with tab_main:
    M_d = mass_slider(st, "Disk kütlesi M_d", DiskParams, "M_d", 1e10, "10¹⁰ M☉", 0.1)
    halo_on = st.toggle("Karanlık madde halesi", value=True)
    inverse = st.toggle("Ters mod: hedef v_c(8.2) → M_vir", value=False, disabled=not halo_on)
    v_target = st.number_input(
        "Hedef v_c(8.2 kpc) [km/s]", 100.0, 400.0, V_TARGET_DEFAULT, 1.0,
        disabled=not (inverse and halo_on),
    )
    M_vir_slider = mass_slider(
        st, "Hale virial kütlesi M_vir", HaloParams, "M_vir", 1e12, "10¹² M☉", 0.05,
        disabled=not halo_on or inverse,
    )
    Q_target = plain_slider(st, "Toomre Q_target", DiskICParams, "Q_target", 0.05)
    lo, hi = param_range(BlackHoleParams, "M_bh")
    log_bh = st.slider(
        "Kara delik log₁₀(M_bh / M☉)", 0.0, float(np.log10(hi)),
        float(np.log10(default(BlackHoleParams, "M_bh"))), 0.01,
        help="Gerçek değer 4.3×10⁶ M☉. Üst sınır bilerek abartılıdır.",
    )
    M_bh = 0.0 if log_bh == 0.0 else 10.0 ** log_bh
    st.caption(f"M_bh = {M_bh:.3g} M☉")

with tab_adv:
    st.markdown("**Disk (Miyamoto–Nagai)**")
    a = plain_slider(st, "a (ölçek uzunluğu)", DiskParams, "a", 0.05, "kpc")
    b = plain_slider(st, "b (ölçek yüksekliği)", DiskParams, "b", 0.01, "kpc")
    st.markdown("**Şişkinlik (Hernquist)**")
    M_b = mass_slider(st, "M_b", BulgeParams, "M_b", 1e10, "10¹⁰ M☉", 0.05)
    c = plain_slider(st, "c", BulgeParams, "c", 0.05, "kpc")
    st.markdown("**Hale (NFW)**")
    r_s = plain_slider(st, "r_s", HaloParams, "r_s", 0.5, "kpc", disabled=not halo_on)
    conc = plain_slider(st, "conc = r_vir / r_s", HaloParams, "conc", 0.5, disabled=not halo_on)
    st.markdown("**Test yıldızı yörüngesi**")
    R0 = plain_slider(st, "Başlangıç yarıçapı R₀", OrbitParams, "R0", 0.1, "kpc")
    v_phi_frac = plain_slider(st, "v_φ / v_c(R₀)", OrbitParams, "v_phi_frac", 0.01)
    v_R0 = plain_slider(st, "v_R", OrbitParams, "v_R0", 1.0, "km/s")
    v_z0 = plain_slider(st, "v_z", OrbitParams, "v_z0", 1.0, "km/s")
    t_orbit = plain_slider(st, "Süre", OrbitParams, "t_orbit", 50.0, "Myr")

with tab_num:
    st.info(
        "Bu parametreler sonucu değiştirmemelidir; "
        "değiştiriyorsa simülasyon yakınsamamıştır."
    )
    dt = plain_slider(st, "dt (zaman adımı)", NumericalParams, "dt", 0.01, "Myr")
    N = int(plain_slider(st, "N (parçacık sayısı)", NumericalParams, "N", 1000.0))
    eps = plain_slider(st, "eps (yumuşatma)", NumericalParams, "eps", 0.005, "kpc")
    theta = plain_slider(st, "theta (Barnes–Hut)", NumericalParams, "theta", 0.05)
    st.caption("Aşama 1'de yalnızca dt kullanılır (yörünge entegrasyonu). N, eps, theta Aşama 2 içindir.")


# ==========================================================================
# Parametreler -> potansiyel
# ==========================================================================

def make_params(M_vir):
    return GalaxyParams(
        disk=DiskParams(M_d=M_d, a=a, b=b),
        bulge=BulgeParams(M_b=M_b, c=c),
        halo=HaloParams(M_vir=M_vir, r_s=r_s, conc=conc, enabled=halo_on),
        bh=BlackHoleParams(M_bh=M_bh),
    )


M_vir = M_vir_slider
if inverse and halo_on:
    sol = solve_mvir_for_vc(make_params(M_vir_slider), v_target, R_SUN)
    if sol is None:
        lo, hi = param_range(HaloParams, "M_vir")
        st.sidebar.error(
            f"{v_target:.0f} km/s için {lo:.0e}–{hi:.0e} M☉ aralığında M_vir çözümü yok. "
            "Kaydırak değeri kullanılıyor."
        )
    else:
        M_vir = sol
        st.sidebar.success(f"Çözülen M_vir = {M_vir:.3e} M☉")

params = make_params(M_vir)
pot = build_potential(params)



# ==========================================================================
# Koşu sonuçları sekmesi (önceden koşulmuş N-cisim koşuları; yalnızca okur)
# ==========================================================================

# Koşu verisi klasörü ve canlı mod parçacık sayısı ortam değişkeniyle değiştirilebilir;
# bulut dağıtımı (deploy/app.py) küçültülmüş veriyi ve daha az parçacığı böyle seçer.
RUNS_DIR = Path(os.environ.get("MW_RUNS_DIR", Path(__file__).resolve().parent / "runs"))


@st.cache_resource(show_spinner=False)
def cached_run(name, stamp):
    """stamp: timeseries.npz değişince önbelleği geçersiz kılmak için mtime."""
    return an.load_run(name, RUNS_DIR)


@st.cache_data(show_spinner="Çubuk ölçümleri hesaplanıyor…")
def cached_bar_history(name, stamp):
    """Çubuk ölçüm tarihçesi. Küçültülmüş (deploy) veride tam çözünürlüklü seri
    bar_history_full.npz olarak gelir; seyreltilmiş snapshot'lardan yeniden hesaplanmaz
    (5000 parçacıkta gürültü tabanı ~2.5 kat yüksek)."""
    run = cached_run(name, stamp)
    full = run.path / "bar_history_full.npz"
    if full.exists():
        return dict(np.load(full))
    return an.bar_history(run)


@st.cache_data(show_spinner=False, max_entries=64)
def cached_snapshot(name, k, stamp):
    return cached_run(name, stamp).snapshot(k)


def _stamp(name):
    return (RUNS_DIR / name / "timeseries.npz").stat().st_mtime


def _step(key, delta, n):
    st.session_state[key] = int(np.clip(st.session_state.get(key, 0) + delta, 0, n - 1))


def render_runs_tab():
    names = an.list_runs(RUNS_DIR)
    if not names:
        st.info("runs/ altında tamamlanmış koşu yok. Koşmak için: python -m src.nbody --name <ad>")
        return
    st.caption("Bu sekme önceden koşulmuş simülasyonları gösterir; kenar çubuğundaki kaydıraklar bunları değiştirmez.")
    name = st.selectbox("Koşu", names, index=names.index("ana") if "ana" in names else 0)
    stamp = _stamp(name)
    run = cached_run(name, stamp)
    st.markdown(f"`{an.run_label(run)}`")
    if "deploy" in run.config:
        d = run.config["deploy"]
        st.caption(f"Küçültülmüş veri: 3B görüntü {d['source_N']} parçacıktan {d['N']}'e seyreltildi, "
                   f"snapshot aralığı {run.config['numerical']['snap_interval']:g} Myr. Çubuk gücü ve enerji "
                   "eğrileri tam çözünürlüklü koşudan.")

    hist = cached_bar_history(name, stamp)
    t_ts, A2 = hist["t"], hist["A2max"]
    present = an.bar_present(hist["A2max"], hist["L_bar"])
    t_bar = an.bar_onset_time(t_ts, present)
    t_E, dE = an.energy_error(run)
    r_com, v_com = an.com_drift(run)
    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("ELN ε", f"{run.config['eln']:.2f}", "< 1.1: çubuk beklenir" if run.config["eln"] < 1.1 else "≥ 1.1",
              delta_color="off")
    m2.metric("Çubuk oluşumu", f"{t_bar:.0f} Myr" if t_bar is not None else "yok",
              f"çubuk var: sürenin %{100 * present.mean():.0f}'i", delta_color="off",
              help=f"Çubuk ölçütü: A2/A0(R)'nin R < {A2_R_MAX:g} kpc kuşaklarındaki maksimumu ≥ {an.BAR_THRESHOLD} "
                   f"ve m=2 fazı en az {an.BAR_MIN_LENGTH:g} kpc boyunca ±{np.degrees(an.BAR_PHASE_TOL):.0f}° sabit. "
                   f"Oluşum zamanı: ölçütün {an.BAR_PERSIST:.0f} Myr kesintisiz sağlandığı ilk an.")
    m3.metric("Maks çubuk gücü", f"{A2.max():.2f}", f"t = {t_ts[A2.argmax()]:.0f} Myr", delta_color="off")
    m4.metric("Maks |ΔE/E|", f"{np.abs(dE).max():.1e}")
    m5.metric("Kütle merkezi kayması", f"{r_com:.3f} kpc", f"maks |v| = {kpcmyr_to_kms(v_com):.1f} km/s",
              delta_color="off")

    key = f"snap_{name}"
    n = run.n_snap
    st.session_state.setdefault(key, 0)
    b1, b2, b3, b4, sl = st.columns([1, 1, 1, 1, 8])
    b1.button("⏮", key=f"{key}_first", on_click=_step, args=(key, -n, n), help="başa")
    b2.button("◀", key=f"{key}_prev", on_click=_step, args=(key, -1, n), help="geri")
    b3.button("▶", key=f"{key}_next", on_click=_step, args=(key, 1, n), help="ileri")
    b4.button("⏭", key=f"{key}_last", on_click=_step, args=(key, n, n), help="sona")
    k = sl.slider("Zaman (snapshot)", 0, n - 1, key=key,
                  format=f"%d  (×{run.config['numerical']['snap_interval']:g} Myr)")
    t = float(run.t_snap[k])
    snap = cached_snapshot(name, k, stamp)
    snap0 = cached_snapshot(name, 0, stamp)
    prof = an.m2_profile(snap, run.mass)
    bm = an.bar_measure(snap, run.mass)

    left, right = st.columns([1.15, 1])
    with left:
        idx = view_subsample(len(run.mass))
        galaxy_3d(snap[idx, :3], UIREV_RUNS, key="runs3d",
                  title=f"{name} — t = {t:.0f} Myr   ({len(idx)}/{len(run.mass)} parçacık çiziliyor)")
        st.caption("Fareyle döndürün, tekerlekle yakınlaştırın. Kamera açısı zaman kaydırağında ve koşu "
                   "değiştirirken korunur.")
    with right:
        overlay = st.checkbox("Diğer koşuların çubuk gücü eğrilerini de göster", value=True)
        fig, ax = plt.subplots(figsize=(6.4, 3.4))
        if overlay:
            for other in names:
                if other != name:
                    o = cached_run(other, _stamp(other))
                    oh = cached_bar_history(other, _stamp(other))
                    ax.plot(oh["t"], an.running_mean(oh["t"], oh["A2max"], 30.0), lw=0.9, alpha=0.6, label=other)
        ax.plot(t_ts, an.running_mean(t_ts, A2, 30.0), color="k", lw=1.6, label=name)
        ax.axvline(t, color="crimson", lw=1)
        ax.axhline(an.BAR_THRESHOLD, color="0.5", ls=":", lw=1)
        ax.set(xlabel="t [Myr]", ylabel="maks A2/A0(R)",
               title=f"Çubuk gücü: A2/A0(R) maksimumu, R < {A2_R_MAX:g} kpc (30 Myr ort.)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)

        fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.8))
        for s_, lab in ((snap0, "t = 0"), (snap, f"t = {t:.0f} Myr")):
            axs[0].semilogy(*an.surface_density_profile(s_, run.mass), label=lab)
        axs[0].set(xlabel="R [kpc]", ylabel="Σ [M☉/kpc²]", title="Yüzey yoğunluğu")
        axs[0].legend(fontsize=7)
        axs[1].plot(t_E, dE, color="0.2", lw=1)
        axs[1].axvline(t, color="crimson", lw=1)
        axs[1].set(xlabel="t [Myr]", title="(E − E₀)/|E₀|")
        for a in axs:
            a.grid(alpha=0.25)
            a.tick_params(labelsize=7)
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
        st.caption(f"R_d = {an.scale_length(snap, run.mass):.2f} kpc (t=0: {an.scale_length(snap0, run.mass):.2f}),  "
                   f"kalınlık medyan|z| = {an.disk_thickness(snap):.3f} kpc (t=0: {an.disk_thickness(snap0):.3f}),  "
                   f"topak kesri = {an.clump_fraction(snap, run.mass):.3f},  "
                   f"çubuk gücü = {bm['A2max']:.2f} @ R = {bm['R_peak']:.2f} kpc, "
                   f"çubuk uzunluğu = {bm['L_bar']:.1f} kpc" if an.bar_present(bm["A2max"], bm["L_bar"])
                   else f"çubuk gücü = {bm['A2max']:.2f} (çubuk ölçütü sağlanmıyor)")

        ok = prof["count"] >= an.BAR_MIN_COUNT
        fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.6))
        axs[0].plot(prof["R_mid"][ok], prof["A2"][ok], "-o", ms=2.5, color="C0")
        axs[0].axhline(an.BAR_THRESHOLD, color="0.5", ls=":", lw=1)
        axs[0].set(xlabel="R [kpc]", title="A2/A0(R)")
        axs[1].axhspan(-np.degrees(an.BAR_PHASE_TOL), np.degrees(an.BAR_PHASE_TOL), color="0.85")
        axs[1].plot(prof["R_mid"][ok], np.degrees(an._dphase(prof["phase"][ok], bm["phase"])), "-o", ms=2.5, color="C1")
        if an.bar_present(bm["A2max"], bm["L_bar"]):
            axs[1].axvspan(0, bm["L_bar"], color="C1", alpha=0.08)
        axs[1].set(xlabel="R [kpc]", title="m=2 fazı − çubuk fazı [°]", ylim=(-95, 95))
        for a in axs:
            a.grid(alpha=0.25)
            a.tick_params(labelsize=7)
            a.title.set_fontsize(9)
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)


# ==========================================================================
# 3B görüntüleme (plotly) — Canlı ve Koşu sonuçları sekmeleri ortak
# ==========================================================================

VIEW3D_BG = "#0b0d14"            # koyu zemin
VIEW3D_HALF = 16.0               # kpc, x ve y eksen yarı genişliği (sabit: sahne zıplamasın)
VIEW3D_Z_HALF = 4.0              # kpc, z eksen yarı genişliği
VIEW3D_R_CMAX = 15.0             # kpc, renk ölçeğinin üst sınırı (yarıçap)
VIEW3D_N = 15000                 # Koşu sonuçları: çizim için seyreltilmiş parçacık sayısı
VIEW3D_SEED = 0                  # seyreltme tohumu: her snapshot'ta aynı parçacıklar
UIREV_RUNS = "kosu-3b"           # sabit uirevision: zaman kaydırağı kamerayı sıfırlamasın
UIREV_LIVE = "canli-3b"


def galaxy_3d_figure(pos, uirevision, title=None, height=600):
    """pos : (N,3) kpc. Yarıçapa göre renkli, küçük yarı saydam noktalar; koyu zemin, minimum eksen.

    uirevision sabit tutulur: aynı değerle yeniden çizildiğinde plotly kullanıcının
    döndürdüğü/yakınlaştırdığı kamera açısını korur.
    """
    R = np.hypot(pos[:, 0], pos[:, 1])
    fig = go.Figure(go.Scatter3d(
        x=pos[:, 0], y=pos[:, 1], z=pos[:, 2], mode="markers", hoverinfo="skip",
        marker=dict(size=1.6, color=R, colorscale="Plasma", cmin=0.0, cmax=VIEW3D_R_CMAX, opacity=0.55,
                    colorbar=dict(title=dict(text="R [kpc]", side="right", font=dict(color="#aab")), thickness=10, len=0.5,
                                  tickfont=dict(color="#aab"), outlinewidth=0)),
    ))
    axis = dict(showgrid=False, zeroline=False, showbackground=False, showticklabels=False,
                showspikes=False, title=dict(text=""), showline=False)
    fig.update_layout(
        uirevision=uirevision,
        paper_bgcolor=VIEW3D_BG, plot_bgcolor=VIEW3D_BG, height=height,
        margin=dict(l=0, r=0, t=30 if title else 0, b=0),
        title=dict(text=title or "", font=dict(color="#dde", size=13), x=0.02),
        scene=dict(
            uirevision=uirevision, bgcolor=VIEW3D_BG,
            xaxis=dict(axis, range=[-VIEW3D_HALF, VIEW3D_HALF]),
            yaxis=dict(axis, range=[-VIEW3D_HALF, VIEW3D_HALF]),
            zaxis=dict(axis, range=[-VIEW3D_Z_HALF, VIEW3D_Z_HALF]),
            aspectmode="manual", aspectratio=dict(x=1, y=1, z=VIEW3D_Z_HALF / VIEW3D_HALF),
            camera=dict(eye=dict(x=0.0, y=-1.25, z=0.85)),
        ),
    )
    return fig


GALAXY3D_DIR = Path(__file__).resolve().parent / "ui" / "galaxy3d"


def _ensure_plotly_js():
    """Bileşenin plotly.js'ini plotly Python paketinden kopyala (internet gerekmez)."""
    import shutil
    import plotly
    dst = GALAXY3D_DIR / "plotly.min.js"
    if not dst.exists():
        shutil.copyfile(Path(plotly.__file__).parent / "package_data" / "plotly.min.js", dst)


_ensure_plotly_js()
_galaxy3d_component = components.declare_component("galaxy3d", path=str(GALAXY3D_DIR))


def _b64_f32(a):
    return base64.b64encode(np.ascontiguousarray(a, dtype=np.float32).tobytes()).decode("ascii")


def galaxy_3d(pos, uirevision, key, title=None, height=600):
    """3B görünümü çiz. st.plotly_chart yerine yerel bileşen: Streamlit 1.65'te st.plotly_chart'ın
    öğe kimliği figür içeriğine bağlı, veri değişince grafik yeniden kurulur ve kamera açısı kaybolur.
    Bu bileşen aynı key ile monte kalır ve Plotly.react ile günceller; uirevision kamerayı korur."""
    spec = json.loads(galaxy_3d_figure(np.zeros((0, 3)), uirevision, title, height).to_json())
    _galaxy3d_component(figure=spec, x=_b64_f32(pos[:, 0]), y=_b64_f32(pos[:, 1]), z=_b64_f32(pos[:, 2]),
                        c=_b64_f32(np.hypot(pos[:, 0], pos[:, 1])), height=height, key=key, default=None)


@st.cache_data(show_spinner=False)
def view_subsample(n_total, n_view=VIEW3D_N, seed=VIEW3D_SEED):
    """Çizim için sabit parçacık alt kümesi (her snapshot'ta aynı indeksler)."""
    if n_total <= n_view:
        return np.arange(n_total)
    return np.sort(np.random.default_rng(seed).choice(n_total, n_view, replace=False))


# ==========================================================================
# Canlı sekmesi: test parçacıkları, analitik potansiyel, gerçek zamanlı leapfrog
# ==========================================================================

LIVE_N = int(os.environ.get("MW_LIVE_N", 5000))  # test parçacığı sayısı
LIVE_SEED = 1
LIVE_FRAME_S = 0.25              # s, kare aralığı (oynarken)
LIVE_R_MIN = 0.5                 # kpc, merkeze çok yakın başlangıç yok
LIVE_ESCAPE_R = 30.0             # kpc, bu yarıçapın dışı "kaçan" sayılır
LIVE_SPEEDS = (1.0, 2.0, 5.0, 10.0)  # Myr / kare


def live_potential():
    """Kenar çubuğundaki fiziksel parametrelerle analitik potansiyel. Kara delik parçacık
    eps'iyle yumuşatılır (merkeze dalan yörüngelerde sayısal patlama olmasın)."""
    return build_potential(params, eps_bh=eps)


def live_reset():
    """Mevcut potansiyelde dairesel yörüngeler: üstel disk R, sech^2 z, v = v_c(R) teğetsel."""
    ic = DiskICParams()
    rng = np.random.default_rng(LIVE_SEED)
    R = np.maximum(sample_radii(rng, LIVE_N, ic.R_d, ic.R_cut), LIVE_R_MIN)
    z = sample_heights(rng, LIVE_N, ic.z_0)
    ph = rng.uniform(0.0, 2.0 * np.pi, LIVE_N)
    v = circular_velocity(live_potential(), R)
    st.session_state.live_pos = np.column_stack((R * np.cos(ph), R * np.sin(ph), z))
    st.session_state.live_vel = np.column_stack((-v * np.sin(ph), v * np.cos(ph), np.zeros(LIVE_N)))
    st.session_state.live_t = 0.0
    st.session_state.live_running = False


def live_set_running(flag):
    st.session_state.live_running = flag


def live_step_once():
    st.session_state.live_step_once = True


def live_advance(n_steps):
    """Mevcut leapfrog (orbits.leapfrog, KDK) ile n_steps adım; potansiyel her karede yeniden kurulur,
    böylece kaydırak değişiklikleri anında etkiler."""
    orb = leapfrog(live_potential(), st.session_state.live_pos, st.session_state.live_vel, dt, n_steps,
                   save_every=n_steps)
    st.session_state.live_pos = orb.pos[-1]
    st.session_state.live_vel = orb.vel[-1]
    st.session_state.live_t += n_steps * dt


def render_live_tab():
    st.caption(f"Test parçacığı rejimi: {LIVE_N} yıldız yalnızca kenar çubuğundaki analitik potansiyelde hareket "
               "eder, birbirlerini çekmezler (öz-yerçekimi yok).")
    if "live_pos" not in st.session_state:
        live_reset()
    c1, c2, c3, c4, c5 = st.columns([1, 1, 1, 1, 3])
    c1.button("▶ Başlat", on_click=live_set_running, args=(True,), disabled=st.session_state.live_running,
              key="live_start")
    c2.button("⏸ Durdur", on_click=live_set_running, args=(False,), disabled=not st.session_state.live_running,
              key="live_stop")
    c3.button("⟲ Sıfırla", on_click=live_reset, key="live_reset",
              help="Parçacıkları mevcut potansiyelde yeniden dairesel yörüngelere koyar.")
    c4.button("⏭ Bir kare", on_click=live_step_once, disabled=st.session_state.live_running, key="live_step")
    speed = c5.select_slider("Hız [Myr / kare]", LIVE_SPEEDS, value=LIVE_SPEEDS[1], key="live_speed")

    @st.fragment(run_every=LIVE_FRAME_S if st.session_state.live_running else None)
    def live_view():
        n_steps = max(1, int(round(speed / dt)))
        if st.session_state.live_running or st.session_state.pop("live_step_once", False):
            live_advance(n_steps)
        pos = st.session_state.live_pos
        r = np.linalg.norm(pos, axis=1)
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Zaman", f"{st.session_state.live_t:.0f} Myr")
        m2.metric(f"Kaçan (r > {LIVE_ESCAPE_R:g} kpc)", f"%{100 * np.mean(r > LIVE_ESCAPE_R):.1f}")
        m3.metric("Medyan yarıçap", f"{np.median(np.hypot(pos[:, 0], pos[:, 1])):.2f} kpc")
        m4.metric("v_c(8.2) şu an", f"{kpcmyr_to_kms(float(circular_velocity(live_potential(), R_SUN))):.0f} km/s")
        galaxy_3d(pos, UIREV_LIVE, key="live3d")

    live_view()
    st.caption(f"Leapfrog, dt = {dt} Myr (Sayısal sekmesi), karede {max(1, int(round(speed / dt)))} adım. "
               "Kaydırakları oynarken değiştirebilirsiniz; yörüngeler yeni potansiyelde devam eder. "
               "Haleyi kapatın: dış disk birkaç yüz Myr içinde dağılır; yıldızlar görüntü kutusunun "
               f"(±{VIEW3D_HALF:g} kpc) dışına çıkar, 'Kaçan' sayacına bakın.")

# ==========================================================================
# Başlık ve sayısal özet
# ==========================================================================

st.title("Samanyolu Dinamiği")
st.markdown("Galaksiyi bir arada tutan şey **görünür madde değil**. Haleyi kapatıp dönme eğrisine bakın.")

tab_model, tab_live, tab_runs = st.tabs(["Analitik model", "Canlı", "Koşu sonuçları"])

with tab_model:

    vc_sun = float(circular_velocity(pot, R_SUN))
    M_dyn = float(dynamical_mass(pot, R_SUN))
    M_true = float(pot.enclosed_mass(R_SUN))

    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("v_c(8.2 kpc)", f"{kpcmyr_to_kms(vc_sun):.1f} km/s")
    c2.metric("Yörünge periyodu", f"{float(orbital_period(pot, R_SUN)):.1f} Myr")
    c3.metric("M_dyn = R·v_c²/G", f"{M_dyn:.3e} M☉")
    c4.metric("M_true (küre içi)", f"{M_true:.3e} M☉", f"M_dyn/M_true = {M_dyn / M_true:.3f}", delta_color="off")
    c5.metric("Kaçış hızı (8.2 kpc)", f"{kpcmyr_to_kms(float(escape_velocity(pot, R_SUN))):.0f} km/s")
    st.caption(
        "M_dyn > M_true çünkü R·v_c²/G küresel simetri varsayar; yassı disk düzlemdeki bir yıldızı "
        "aynı kütleli bir küreden daha güçlü çeker — küresel bir dağılımda kabuk teoremi gereği oran tam 1.00 olurdu."
    )


    @st.cache_data(show_spinner="N-cisim diski üretiliyor…")
    def drift_summary(bulge, halo, bh, Q_target, num):
        """Varsayılan N-cisim diski (M_d = DiskICParams.M_d) için 8.2 kpc sürüklenme tanısı, km/s."""
        gp = GalaxyParams(bulge=bulge, halo=halo, bh=bh)
        ic = DiskICParams(Q_target=Q_target)
        disk = make_disk(gp, ic, num, verbose=False)
        return asymmetric_drift_check(gp, ic, disk, R_SUN, DRIFT_HALF_WIDTH), ic.M_d


    drift, M_d_nbody = drift_summary(params.bulge, params.halo, params.bh, Q_target,
                                     NumericalParams(N=N, eps=eps, dt=dt, theta=theta))
    d1, d2, d3, d4 = st.columns(4)
    d1.metric("v_c(8.2), N-cisim parçacıkları", f"{drift['v_c']:.1f} km/s")
    d2.metric("⟨v_φ⟩(8.2), N-cisim parçacıkları", f"{drift['vphi_mean']:.1f} km/s")
    d3.metric("Sürüklenme v_c − ⟨v_φ⟩, ölçülen", f"{drift['measured']:.2f} ± {drift['measured_err']:.2f} km/s")
    d4.metric("Sürüklenme, Jeans denklemi", f"{drift['predicted']:.2f} km/s",
              f"ölçülen/öngörülen = {drift['measured'] / drift['predicted']:.3f}", delta_color="off")
    st.caption(
        "Gözlem soğuk izleyicilerden (gaz, genç yıldızlar, maserler) yapılır; onların ⟨v_φ⟩'si v_c'ye eşittir, "
        "bu yüzden 229 km/s dairesel hızla karşılaştırılır — sıcak disk yıldızları ise sürüklenme kadar yavaş döner. "
        f"(N-cisim diski: M_d = {M_d_nbody:.3e} M☉, Q_target ve hale kaydıraklardan, R = 8.2 ± {DRIFT_HALF_WIDTH} kpc.)"
    )


    # ==========================================================================
    # Panel 1 & 2 — dönme eğrileri
    # ==========================================================================

    def plot_curves(ax, R):
        comps = circular_velocity_components(pot, R)
        for name, v in comps.items():
            ax.plot(R, kpcmyr_to_kms(v), color=COLORS[name], lw=1.4, ls="--", label=LABELS[name])
        ax.plot(R, kpcmyr_to_kms(circular_velocity(pot, R)), color=COLORS["total"], lw=2.4, label=LABELS["total"])
        return comps


    col_a, col_b = st.columns(2)

    with col_a:
        st.subheader("Dönme eğrisi")
        R = np.linspace(*PLOT_R_LIN, 600)
        fig, ax = plt.subplots(figsize=(6.4, 4.3))
        plot_curves(ax, R)
        o = OBS_EILERS2019
        R_obs = np.linspace(o["R_min"], o["R_max"], o["n"])
        ax.errorbar(R_obs, o["v0"] + o["slope"] * (R_obs - o["R0"]), yerr=o["err"],
                    fmt="o", ms=4, color="#e7298a", label="Eilers+2019 uyum eğrisi")
        ax.axvline(R_SUN, color="0.6", lw=0.8, ls=":")
        ax.set(xlabel="R [kpc]", ylabel="v_c [km/s]", xlim=PLOT_R_LIN, ylim=(0, 320))
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(alpha=0.25)
        st.pyplot(fig)
        plt.close(fig)

    with col_b:
        st.subheader("Dört rejim (logaritmik yarıçap)")
        R = np.logspace(np.log10(PLOT_R_LOG[0]), np.log10(PLOT_R_LOG[1]), 800)
        fig, ax = plt.subplots(figsize=(6.4, 4.3))
        comps = plot_curves(ax, R)
        # Her yarıçapta v_c^2'ye en çok katkı veren bileşenin bölgesini boya
        names = list(comps)
        dom = np.argmax(np.vstack([comps[n] for n in names]), axis=0)
        edges = np.flatnonzero(np.diff(dom)) + 1
        for s, e in zip(np.r_[0, edges], np.r_[edges, len(R) - 1]):
            ax.axvspan(R[s], R[e], color=COLORS[names[dom[s]]], alpha=0.08, lw=0)
        ax.set(xscale="log", yscale="log", xlabel="R [kpc]", ylabel="v_c [km/s]",
               xlim=PLOT_R_LOG, ylim=(1, 1000))
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(alpha=0.25, which="both")
        st.pyplot(fig)
        plt.close(fig)
        st.caption("Arka plan rengi, o yarıçapta dönme hızına en çok katkı veren bileşeni gösterir.")


    # ==========================================================================
    # Panel 4 — Q(R) ve test yıldızı yörüngesi
    # ==========================================================================

    col_c, col_d = st.columns(2)

    with col_c:
        st.subheader("Toomre Q(R)")
        st.markdown("*varsayılan Sigma ve sigma_R ile*")
        ic = DiskICParams(Q_target=Q_target)
        Sigma_fn = lambda r: exponential_surface_density(r, M_d, ic.R_d, ic.R_cut)
        R = np.linspace(0.5, ic.R_cut * 0.999, 300)
        sig = sigma_R_profile(pot, R, Sigma_fn, Q_target, Q_REF_RADIUS, ic.R_d)
        Q = toomre_q(pot, R, sig, Sigma_fn(R))
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
        ax.plot(R, Q, color=COLORS["disk"], lw=2)
        ax.axhline(1.0, color="crimson", lw=1, ls="--", label="Q = 1 (kararlılık sınırı)")
        ax.fill_between(R, 0, 1, color="crimson", alpha=0.06)
        ax.axvline(Q_REF_RADIUS, color="0.6", lw=0.8, ls=":")
        ax.set(xlabel="R [kpc]", ylabel="Q", ylim=(0, max(3.0, 1.1 * np.nanmax(Q[R > 2]))))
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
        st.pyplot(fig)
        plt.close(fig)
        st.caption(
            f"Üstel disk (R_d = {ic.R_d} kpc), σ_R ∝ exp(−R/2R_d), Q({Q_REF_RADIUS} kpc) = Q_target. "
            f"σ_R({Q_REF_RADIUS} kpc) = {kpcmyr_to_kms(float(sigma_R_profile(pot, Q_REF_RADIUS, Sigma_fn, Q_target, Q_REF_RADIUS, ic.R_d))):.1f} km/s."
        )


    @st.cache_data(show_spinner="Yörünge entegre ediliyor…")
    def run_orbit(params, R0, v_phi_frac, v_R0, v_z0, t_orbit, dt):
        p = build_potential(params)
        pos, vel = initial_conditions(p, R0, v_phi_frac, v_R0, v_z0)
        n_steps = int(round(t_orbit / dt))
        orb = leapfrog(p, pos, vel, dt, n_steps, save_every=max(1, n_steps // MAX_ORBIT_SAMPLES))
        E = energy(p, orb.pos[:, 0], orb.vel[:, 0])
        Lz = angular_momentum_z(orb.pos[:, 0], orb.vel[:, 0])
        return orb.t, orb.pos[:, 0], E, Lz


    with col_d:
        st.subheader("Test yıldızı yörüngesi (leapfrog)")
        t, xyz, E, Lz = run_orbit(params, R0, v_phi_frac, v_R0, v_z0, t_orbit, dt)
        fig, axs = plt.subplots(1, 2, figsize=(6.4, 3.2))
        axs[0].plot(xyz[:, 0], xyz[:, 1], lw=0.5, color=COLORS["halo"])
        axs[0].plot(0, 0, "k+")
        axs[0].set(xlabel="x [kpc]", ylabel="y [kpc]", aspect="equal", title="Disk düzlemi")
        axs[1].plot(np.hypot(xyz[:, 0], xyz[:, 1]), xyz[:, 2], lw=0.5, color=COLORS["halo"])
        axs[1].set(xlabel="R [kpc]", ylabel="z [kpc]", title="Meridyen düzlemi")
        for x in axs:
            x.grid(alpha=0.25)
            x.tick_params(labelsize=8)
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.4, 2.2))
        ax.plot(t, (E - E[0]) / abs(E[0]), lw=0.8, color="0.2", label="ΔE/|E|")
        ax.plot(t, (Lz - Lz[0]) / abs(Lz[0]) if Lz[0] != 0 else np.zeros_like(t),
                lw=0.8, color=COLORS["bulge"], label="ΔL_z/|L_z|")
        ax.set(xlabel="t [Myr]", ylabel="bağıl sapma")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
        ax.tick_params(labelsize=8)
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
        st.caption(
            f"{len(t) - 1} kayıt, dt = {dt} Myr. max|ΔE/E| = {np.max(np.abs(E - E[0]) / abs(E[0])):.1e} — "
            "leapfrog simplektik olduğu için hata salınır ama büyümez."
        )


with tab_live:
    render_live_tab()


with tab_runs:
    render_runs_tab()
