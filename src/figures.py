"""Figures for the README.

needles.png  every tile's giant-dune crest drawn as a needle on the map, coloured by
             which climate's predicted crest it sits closer to
dials.png    one dial per dune field: observed crest directions against the modern
             and LGM predicted crests
tests.png    every test's result against what pure-modern and perfect-memory crests
             would give
"""
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from addendum import base_tiles
from analysis import load, pred
from flux import axial_diff, signed_axial

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
RES = ROOT / "results"

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
SURFACE, GRID = "#fcfcfb", "#e4e3df"
MODERN, LGM, NEUTRAL = "#2a78d6", "#eb6834", "#b9b8b2"
LAND, SEA = "#f3f1ea", "#fcfcfb"

mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "savefig.dpi": 200,
})

REGIONS = [
    ("North Africa and Arabia", (-18, 60, 12, 35)),
    ("Central and East Asia", (55, 108, 25, 45)),
    ("Simpson Desert (Munga-Thirri)", (134.5, 140.5, -27.8, -22.2)),
    ("Namib", (13.5, 17.5, -28, -22)),
]


def needle(ax, lon, lat, az, length, **kw):
    """Line through (lon, lat) along compass azimuth az, length in degrees of latitude."""
    a = np.radians(az)
    dx = 0.5 * length * np.sin(a) / np.cos(np.radians(lat))
    dy = 0.5 * length * np.cos(a)
    ax.plot([lon - dx, lon + dx], [lat - dy, lat + dy], transform=ccrs.PlateCarree(),
            solid_capstyle="round", **kw)


def classify(d, past="LGM_ens", margin=5.0):
    m_mod = axial_diff(d.large_crest, pred(d, "mod", "type"))
    m_past = axial_diff(d.large_crest, pred(d, past, "type"))
    return np.where(m_mod + margin < m_past, "modern", np.where(m_past + margin < m_mod, "lgm", "tie"))


def fig_needles(d):
    cls = classify(d)
    colour = {"modern": MODERN, "lgm": LGM, "tie": NEUTRAL}
    fig = plt.figure(figsize=(12, 8.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], width_ratios=[1.6, 1, 0.75],
                          hspace=0.16, wspace=0.08, top=0.86, bottom=0.02, left=0.02, right=0.98)
    slots = [gs[0, :], gs[1, 0], gs[1, 1], gs[1, 2]]
    for (name, ext), slot in zip(REGIONS, slots):
        ax = fig.add_subplot(slot, projection=ccrs.PlateCarree())
        ax.set_extent(ext, ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor=LAND, edgecolor="none")
        ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor=SEA, edgecolor="none")
        ax.add_feature(cfeature.COASTLINE.with_scale("50m"), edgecolor=MUTED, linewidth=0.4)
        ax.add_feature(cfeature.BORDERS.with_scale("50m"), edgecolor=GRID, linewidth=0.4)
        span = ext[3] - ext[2]
        length = 0.28 if span < 12 else 0.32
        inside = (d.lon.between(ext[0], ext[1])) & (d.lat.between(ext[2], ext[3]))
        for r, c in zip(d[inside].itertuples(), cls[inside.values]):
            needle(ax, r.lon, r.lat, r.large_crest, length, color=colour[c],
                   linewidth=1.3 if c != "tie" else 0.9, alpha=0.95 if c != "tie" else 0.8)
        ax.set_title(name, loc="left", fontsize=10, color=INK, pad=4)
        ax.set_anchor("NW")
        for s in ax.spines.values():
            s.set_edgecolor(GRID)
    n = pd.Series(cls).value_counts()
    handles = [Line2D([], [], color=MODERN, lw=2, label=f"closer to today's winds ({n.get('modern', 0)})"),
               Line2D([], [], color=LGM, lw=2, label=f"closer to ice-age winds ({n.get('lgm', 0)})"),
               Line2D([], [], color=NEUTRAL, lw=2, label=f"within 5\u00b0 of both ({n.get('tie', 0)})")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.02, 0.925), frameon=False,
               fontsize=9, ncol=3, handlelength=1.6, columnspacing=2.0)
    fig.suptitle("Which wind built the giant dunes?", x=0.02, ha="left",
                 fontsize=15, color=INK, y=0.995)
    fig.text(0.02, 0.945, "One needle per 31 km tile, along the measured crest of the largest dunes. "
             "Colour: is the crest closer to the direction predicted from today's winds (ERA5) or from "
             "Last Glacial Maximum winds (10-model ensemble)?", fontsize=9, color=INK2)
    fig.savefig(FIG / "needles.png", bbox_inches="tight")
    plt.close(fig)


def axial_rose(ax, angles, color, bins=36):
    a = np.radians(np.r_[angles, np.asarray(angles) + 180.0])
    h, e = np.histogram(a, bins=bins * 2, range=(0, 2 * np.pi))
    h = h / h.max()
    ax.bar(e[:-1] + np.pi / bins / 2, h, width=2 * np.pi / (bins * 2) * 0.9, color=color,
           edgecolor="none", alpha=0.9)


def circ_axial_median(x):
    z = np.nanmean(np.exp(2j * np.radians(x)))
    return (np.degrees(np.angle(z)) / 2.0) % 180.0


def fig_dials(d):
    counts = d.field.value_counts()
    fields = [f for f in counts.index if counts[f] >= 12]
    fields = sorted(fields, key=lambda f: d[d.field == f].lon.median())
    ncol = 6
    nrow = int(np.ceil(len(fields) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(12, 2.3 * nrow + 0.8),
                             subplot_kw={"projection": "polar"})
    for ax in axes.ravel():
        ax.set_visible(False)
    for ax, f in zip(axes.ravel(), fields):
        s = d[d.field == f]
        ax.set_visible(True)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        axial_rose(ax, s.large_crest.values, "#c9c7bf")
        for col, colour in (("mod", MODERN), ("LGM_ens", LGM)):
            m = circ_axial_median(pred(s, col, "type"))
            for t in (m, m + 180):
                ax.plot([np.radians(t)] * 2, [0, 1.12], color=colour, lw=2, solid_capstyle="round")
        ax.set_ylim(0, 1.15)
        ax.set_yticks([])
        ax.set_xticks(np.radians([0, 90, 180, 270]))
        ax.set_xticklabels(["N", "E", "S", "W"], fontsize=7, color=MUTED)
        ax.grid(color=GRID, linewidth=0.5)
        ax.spines["polar"].set_edgecolor(GRID)
        ax.set_title(f"{f}\n{len(s)} tiles", fontsize=8, color=INK, pad=6)
    handles = [mpl.patches.Patch(color="#c9c7bf", label="measured giant-dune crests"),
               Line2D([], [], color=MODERN, lw=2, label="crest predicted from today's winds"),
               Line2D([], [], color=LGM, lw=2, label="crest predicted from LGM winds")]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG / "dials.png", bbox_inches="tight")
    plt.close(fig)


def fig_tests():
    t = pd.read_csv(RES / "addendum.csv")
    show = ["Primary: all tiles", "Primary: informative", "B: 10 models agree",
            "A: daily winds, 4 m/s", "A: daily winds, 4 m/s, informative",
            "A: daily winds, 4 m/s, AWIESM11LR", "A: daily winds, 4 m/s, MPIESM12LR",
            "B: both daily models agree"]
    names = {"Primary: all tiles": "Monthly winds, 10 models: all tiles",
             "Primary: informative": "Monthly winds: tiles where climates differ",
             "B: 10 models agree": "Monthly winds: tiles where models agree",
             "A: daily winds, 4 m/s": "Daily winds, 2 models: all tiles",
             "A: daily winds, 4 m/s, informative": "Daily winds: tiles where climates differ",
             "A: daily winds, 4 m/s, AWIESM11LR": "Daily winds, AWI-ESM only",
             "A: daily winds, 4 m/s, MPIESM12LR": "Daily winds, MPI-ESM only",
             "B: both daily models agree": "Daily winds: tiles where both models agree"}
    t = t.set_index("test").loc[show].reset_index()
    t["test"] = t.test.map(names)
    y = np.arange(len(t))[::-1]
    fig, ax = plt.subplots(figsize=(9, 0.55 * len(t) + 1.6))
    ax.axvline(0, color=GRID, lw=1, zorder=0)
    for yi, r in zip(y, t.itertuples()):
        ax.plot([r.lo, r.hi], [yi, yi], color=INK2, lw=1.5, solid_capstyle="round")
        ax.plot(r.expected_if_modern, yi, marker="|", ms=14, mew=2.5, color=MODERN)
        ax.plot(r.expected_if_memory, yi, marker="|", ms=14, mew=2.5, color=LGM)
        ax.plot(r.gain, yi, "o", ms=8, color=INK, mec=SURFACE, mew=1.5, zorder=3)
        ax.text(1.01, yi, f"{r.n} tiles, {r.fields} fields", transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=7.5, color=MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels(t.test, fontsize=8.5, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.7, len(t) - 0.3)
    lim = max(20, np.nanmax(np.abs(t[["lo", "hi"]].values)) + 2)
    ax.set_xlim(-lim, lim)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("degrees by which LGM winds fit giant crests better than today's winds")
    handles = [Line2D([], [], marker="o", ls="", color=INK, ms=7, label="observed (95% field-bootstrap range)"),
               Line2D([], [], marker="|", ls="", color=MODERN, ms=12, mew=2.5, label="if crests follow today's winds"),
               Line2D([], [], marker="|", ls="", color=LGM, ms=12, mew=2.5, label="if crests remember the LGM")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=3, frameon=False,
              fontsize=8, handletextpad=0.4, columnspacing=1.6)
    fig.tight_layout()
    fig.savefig(FIG / "tests.png", bbox_inches="tight")
    plt.close(fig)


def fig_models(d):
    """Do the models agree on which way the LGM turned the crests?"""
    from addendum import add_daily, signed_rot
    dd, _ = add_daily(load())
    b = base_tiles(dd)
    a = signed_rot(b, "LGMd4_AWIESM11LR")
    m = signed_rot(b, "LGMd4_MPIESM12LR")
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    lim = 90
    ax.axhspan(0, lim, xmin=0.5, xmax=1, color="#f0efec", zorder=0)
    ax.axhspan(-lim, 0, xmin=0, xmax=0.5, color="#f0efec", zorder=0)
    ax.axhline(0, color=GRID, lw=1)
    ax.axvline(0, color=GRID, lw=1)
    ax.scatter(a, m, s=9, color=INK2, alpha=0.35, edgecolor="none")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xticks([-90, -45, 0, 45, 90])
    ax.set_yticks([-90, -45, 0, 45, 90])
    ax.set_xlabel("AWI-ESM: LGM turn of the predicted crest (\u00b0)")
    ax.set_ylabel("MPI-ESM: LGM turn of the predicted crest (\u00b0)")
    r = np.corrcoef(a, m)[0, 1]
    same = np.mean(np.sign(a) == np.sign(m))
    ax.set_title(f"Two models, same tiles: same direction at {same:.0%} of tiles, r = {r:.2f}",
                 fontsize=9, color=INK, loc="left")
    ax.text(85, 85, "models agree", ha="right", va="top", fontsize=8, color=MUTED)
    ax.text(-85, -85, "models agree", ha="left", va="bottom", fontsize=8, color=MUTED)
    fig.tight_layout()
    fig.savefig(FIG / "models.png", bbox_inches="tight")
    plt.close(fig)


def fig_turnover(d):
    """Misfit to today's winds against (a) cross-section turnover time and
    (b) the time to turn the whole dune, pooled and within fields."""
    from reorient import t_rot
    from turnover import prepare
    lin = d[d.large_type.eq("l")]
    panels = [
        (d, "a  Time to rebuild the cross-section, all giant dunes",
         prepare(d, d.large_height.values, d.large_wl.values),
         axial_diff(d.large_crest.values, pred(d, "mod", "type")), "crest mode set by dune type"),
        (lin, "b  Time to turn the whole dune, linear dunes", np.log10(t_rot(lin)),
         axial_diff(lin.large_crest.values, lin.mod_f.values), "elongation mode"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, (dd, title, logT, mis, note) in zip(axes, panels):
        ok = np.isfinite(logT)
        logT, mis, fields = logT[ok], np.asarray(mis)[ok], dd.field.values[ok]
        within = mis - pd.Series(mis).groupby(fields).transform("mean").values + mis.mean()
        ax.scatter(10 ** logT, mis, s=6, color="#c9c7bf", edgecolor="none", alpha=0.6, zorder=1)
        edges = np.quantile(logT, np.linspace(0, 1, 9))
        mid = 10 ** ((edges[:-1] + edges[1:]) / 2)
        for vals, colour, lab in ((mis, INK2, "all tiles, binned"),
                                  (within, MODERN, "within each sand sea, binned")):
            m = [vals[(logT >= a) & (logT <= b)].mean() for a, b in zip(edges[:-1], edges[1:])]
            ax.plot(mid, m, "-o", color=colour, lw=2, ms=6, mec=SURFACE, mew=1.2, label=lab, zorder=3)
        ax.axvline(21000, color=LGM, lw=2, zorder=2)
        share = np.mean(10 ** logT > 21000)
        ax.text(21000 * 1.1, 86, "LGM, 21,000" + chr(10) + f"years ago" + chr(10) + f"{share:.0%} of tiles" + chr(10) + "beyond",
                color=LGM, fontsize=8, va="top")
        # chance: crests shuffled between tiles, scored the same way
        rng = np.random.default_rng(0)
        obs = dd.large_crest.values[ok]
        mt = pred(dd, "mod", "type")[ok] if "type" in note else dd.mod_f.values[ok]
        score = lambda c: axial_diff(c, mt)
        chance = np.mean([score(rng.permutation(obs)).mean() for _ in range(200)])
        ax.axhline(chance, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=0)
        ax.text(110, chance + 1.2, f"chance ({chance:.0f}°, crests shuffled)", color=MUTED, fontsize=7.5)
        ax.set_xscale("log")
        ax.set_xlim(100, 3e5)
        ax.set_ylim(0, 90)
        ax.set_yticks([0, 15, 30, 45, 60, 75, 90])
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        ax.text(0.01, 0.01, note, transform=ax.transAxes, fontsize=7.5, color=MUTED)
        ax.set_xlabel("years")
    axes[0].set_ylabel("giant crest vs today's winds (°)")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "turnover.png", bbox_inches="tight")
    plt.close(fig)


def pick_showcase(d):
    """A clear giant-dune tile where the models predict a sizeable LGM rotation
    (chosen from predictions and crest clarity only, for illustration)."""
    rot = np.abs(signed_axial(pred(d, "LGM_ens", "type"), pred(d, "mod", "type")))
    cand = d[(rot >= 15) & (d.large_R >= 0.8)]
    return cand.sort_values("large_contrast", ascending=False).iloc[0]


def fig_tile(d):
    import xarray as xr
    from matplotlib.colors import LightSource
    import orient
    from flux import flux_rose
    from predict import MODELS, model_delta

    r = pick_showcase(d)
    z = orient.read_dem(r.lat, r.lon, half=0.07)          # ~15 km window, streamed
    dx = orient.ARCSEC_M * np.cos(np.radians(r.lat))
    ls = LightSource(azdeg=315, altdeg=35)
    shade = ls.hillshade(z, vert_exag=3, dx=dx, dy=orient.ARCSEC_M)

    era = np.load(ROOT / "data" / "processed" / "era5_winds.npz")
    u, v = era["u"][int(r.tile)].astype(float), era["v"][int(r.tile)].astype(float)
    month = pd.to_datetime(era["hours"]).month.values - 1
    dus, dvs = zip(*[model_delta(m, "lgm", np.array([r.lat]), np.array([r.lon])) for m in MODELS])
    du, dv = np.mean(dus, 0)[0], np.mean(dvs, 0)[0]
    rose_mod = flux_rose(u[None], v[None])[0]
    rose_lgm = flux_rose((u + du[month])[None], (v + dv[month])[None])[0]

    fig = plt.figure(figsize=(11, 5.4))
    ax = fig.add_axes([0.02, 0.06, 0.5, 0.82])
    ny, nx = z.shape
    ext_km = (0, nx * dx / 1000, 0, ny * orient.ARCSEC_M / 1000)
    ax.imshow(shade, cmap="gray", extent=ext_km, vmin=0.15, vmax=1.0)
    cx, cy = ext_km[1] / 2, ext_km[3] / 2
    L = 0.42 * min(ext_km[1], ext_km[3])
    import matplotlib.patheffects as pe
    halo = [pe.Stroke(linewidth=5.2, foreground=INK), pe.Normal()]
    for az, colour, lw, ls_ in ((r.large_crest, "#ffffff", 3.2, "-"),
                                (pred(d.loc[[r.name]], "mod", "type")[0], MODERN, 2.4, "-"),
                                (pred(d.loc[[r.name]], "LGM_ens", "type")[0], LGM, 2.4, "-")):
        a = np.radians(az)
        ax.plot([cx - L * np.sin(a), cx + L * np.sin(a)], [cy - L * np.cos(a), cy + L * np.cos(a)],
                color=colour, lw=lw, ls=ls_, solid_capstyle="round",
                path_effects=halo if colour == "#ffffff" else None)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(0, 1.02, f"{r.field}, {abs(r.lat):.2f}\u00b0{'N' if r.lat >= 0 else 'S'} "
            f"{abs(r.lon):.2f}\u00b0{'E' if r.lon >= 0 else 'W'}  \u00b7  {ext_km[1]:.0f} km across",
            transform=ax.transAxes, fontsize=9, color=INK2)

    coarse = lambda x: x.reshape(36, 5).sum(1)          # 10 degree bins for display
    rose_mod, rose_lgm = coarse(rose_mod), coarse(rose_lgm)
    th = np.radians((np.arange(36) + 0.5) * 10)
    top = max(rose_mod.max(), rose_lgm.max())
    for i, (rose, colour, title) in enumerate(((rose_mod, MODERN, "Sand flux today (ERA5)"),
                                               (rose_lgm, LGM, "Sand flux at the LGM (10 models)"))):
        pax = fig.add_axes([0.57 + i * 0.235, 0.22, 0.17, 0.5], projection="polar")
        pax.set_theta_zero_location("N")
        pax.set_theta_direction(-1)
        pax.bar(th, rose / top, width=np.radians(10) * 0.85, color=colour, edgecolor="none")
        pax.set_yticks([])
        pax.set_xticks(np.radians([0, 90, 180, 270]))
        pax.set_xticklabels(["N", "E", "S", "W"], fontsize=7, color=MUTED)
        pax.grid(color=GRID, lw=0.5)
        pax.spines["polar"].set_edgecolor(GRID)
        pax.set_title(title, fontsize=9, color=INK, pad=10)
    handles = [Line2D([], [], color="#ffffff", lw=3, label="measured crest of the giant dunes",
                      path_effects=[pe.Stroke(linewidth=5, foreground=INK), pe.Normal()]),
               Line2D([], [], color=MODERN, lw=2.4, label="crest predicted from today's winds"),
               Line2D([], [], color=LGM, lw=2.4, label="crest predicted from LGM winds")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.56, 0.02), frameon=False, fontsize=8.5)
    fig.text(0.56, 0.86, "Roses show where sand is blown to, by direction.", fontsize=8.5, color=INK2)
    fig.savefig(FIG / "tile.png", bbox_inches="tight")
    plt.close(fig)
    return r


def main():
    FIG.mkdir(exist_ok=True)
    d = base_tiles(load())
    fig_needles(d)
    fig_dials(d)
    fig_tests()
    fig_models(d)
    fig_turnover(d)
    fig_tile(d)
    print("figures written")


if __name__ == "__main__":
    main()
