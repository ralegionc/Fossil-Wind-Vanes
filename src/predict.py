"""Predicted crest orientations at every tile for the modern, mid-Holocene and LGM climates.

Modern: flux roses from the sampled ERA5 hours.
Past: the delta method. Each model's monthly mean wind change (past run minus its
pre-industrial control) is added to every ERA5 hour of the same calendar month,
so model biases in the mean climate cancel and only the simulated change is used.
The ensemble prediction uses the multi-model mean change; each model is also kept.

Uncertainty in the modern prediction from hour sampling is estimated by a
bootstrap over sampled days.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from flux import axial_diff, flux_rose, orientation_bi, orientation_f, resultant

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
MODELS = ["AWI-ESM-1-1-LR", "INM-CM4-8", "MIROC-ES2L", "MPI-ESM1-2-LR",
          "CNRM-CM5", "GISS-E2-R", "IPSL-CM5A-LR", "MIROC-ESM", "MPI-ESM-P", "MRI-CGCM3"]
EPOCHS = {"lgm": "LGM", "midHolocene": "MH"}
N_BOOT = 100


def control_for(model, epoch):
    return "piControl_p141" if (model == "GISS-E2-R" and epoch == "lgm") else "piControl"


def at_tiles(da, lat, lon):
    """Bilinear values of a (month, lat, lon) field at tile points, wrapping longitude."""
    da = da.assign_coords(lon=(da.lon % 360.0)).sortby("lon").sortby("lat")
    pad = xr.concat([da.isel(lon=[-1]).assign_coords(lon=da.lon[-1:] - 360.0), da,
                     da.isel(lon=[0]).assign_coords(lon=da.lon[:1] + 360.0)], dim="lon")
    pts = pad.interp(lat=xr.DataArray(lat, dims="tile"), lon=xr.DataArray(lon % 360.0, dims="tile"))
    return pts.transpose("tile", "month").values          # (n_tiles, 12)


def model_delta(model, epoch, lat, lon):
    out = []
    for var in ("uas", "vas"):
        a = xr.open_dataset(PROC / "pmip" / f"{model}_{epoch}_{var}.nc")[var]
        b = xr.open_dataset(PROC / "pmip" / f"{model}_{control_for(model, epoch)}_{var}.nc")[var]
        out.append(at_tiles(a - b, lat, lon))
    return out  # du, dv, each (n_tiles, 12)


def predictions(u, v, tag):
    rose = flux_rose(u, v)
    rdd, rdp, dp = resultant(rose)
    return {f"{tag}_bi": orientation_bi(rose), f"{tag}_f": orientation_f(rose),
            f"{tag}_rdd": rdd, f"{tag}_rdp": rdp, f"{tag}_dp": dp}


def axial_sd(x):
    """Circular standard deviation (deg) of axial angles along the last axis."""
    z = np.nanmean(np.exp(2j * np.radians(x)), axis=-1)
    return np.degrees(np.sqrt(-2.0 * np.log(np.abs(z) + 1e-12)) / 2.0)


def main():
    tiles = pd.read_csv(PROC / "tiles.csv")
    era = np.load(PROC / "era5_winds.npz")
    u, v = era["u"].astype(np.float64), era["v"].astype(np.float64)
    hours = pd.to_datetime(era["hours"])
    month = hours.month.values - 1
    lat, lon = tiles.lat.values, tiles.lon.values

    out = pd.DataFrame({"tile": tiles.tile})
    for k, val in predictions(u, v, "mod").items():
        out[k] = val

    # bootstrap over days for the modern prediction
    rng = np.random.default_rng(1)
    day = hours.normalize().values
    days, inv = np.unique(day, return_inverse=True)
    members = [np.flatnonzero(inv == i) for i in range(len(days))]
    boot_bi, boot_f = [], []
    for _ in range(N_BOOT):
        pick = np.concatenate([members[i] for i in rng.integers(0, len(days), len(days))])
        rose = flux_rose(u[:, pick], v[:, pick])
        boot_bi.append(orientation_bi(rose))
        boot_f.append(orientation_f(rose))
    out["mod_bi_sd"] = axial_sd(np.stack(boot_bi, -1))
    out["mod_f_sd"] = axial_sd(np.stack(boot_f, -1))

    for epoch, tag in EPOCHS.items():
        dus, dvs = [], []
        for m in MODELS:
            try:
                du, dv = model_delta(m, epoch, lat, lon)
            except FileNotFoundError:
                print("missing", m, epoch)
                continue
            dus.append(du)
            dvs.append(dv)
            short = m.replace("-", "")
            for k, val in predictions(u + du[:, month], v + dv[:, month], f"{tag}_{short}").items():
                if k.endswith(("_bi", "_f")):
                    out[k] = val
            out[f"{tag}_{short}_du"] = du.mean(1)
            out[f"{tag}_{short}_dv"] = dv.mean(1)
        du, dv = np.mean(dus, 0), np.mean(dvs, 0)
        for k, val in predictions(u + du[:, month], v + dv[:, month], f"{tag}_ens").items():
            out[k] = val
        out[f"{tag}_ens_du"], out[f"{tag}_ens_dv"] = du.mean(1), dv.mean(1)
        out[f"{tag}_nmodels"] = len(dus)

    out.to_csv(PROC / "predictions.csv", index=False)
    print(out.shape)
    for tag in EPOCHS.values():
        print(tag, "median rotation of BI crest vs modern:",
              np.nanmedian(axial_diff(out[f"{tag}_ens_bi"], out["mod_bi"])).round(1),
              " F:", np.nanmedian(axial_diff(out[f"{tag}_ens_f"], out["mod_f"])).round(1))
    print("modern BI bootstrap sd median", out.mod_bi_sd.median().round(1), " F", out.mod_f_sd.median().round(1))


if __name__ == "__main__":
    main()
