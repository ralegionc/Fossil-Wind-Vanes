"""Load the Gunn et al. (2022) per-tile dune geometry table.

Each row is one ~31 km ERA5 grid tile inside a dune field, with a 'large' and a
'small' dune type and wavelength. Where only one dune scale was found the two
are duplicates.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "gunn2022" / "source_data.xlsx"
OUT = ROOT / "data" / "processed" / "tiles.csv"

COLS = {
    "Dune Field": "field",
    "Latitude": "lat",
    "Longitude": "lon",
    "Large Type": "large_type",
    "Large Wavelength (m)": "large_wl",
    "Large Width (m)": "large_width",
    "Large Height (m)": "large_height",
    "Small Type": "small_type",
    "Small Wavelength (m)": "small_wl",
    "Small Width (m)": "small_width",
    "Small Height (m)": "small_height",
    "Mean Resultant Flux Magnitude (m^2/yr)": "rdp",
    "Mean of Flux Magnitudes (m^2/yr)": "dp",
}


def load_tiles() -> pd.DataFrame:
    df = pd.read_excel(RAW, sheet_name="Dune Geometries", header=2)
    df = df.rename(columns=COLS)[list(COLS.values())].dropna(subset=["lat"])
    df["field"] = df["field"].str.replace("\\&", "&", regex=False)
    for c in df.columns:
        if c not in ("field", "large_type", "small_type"):
            df[c] = df[c].astype(float)
    df["lon"] = (df["lon"] + 180.0) % 360.0 - 180.0  # some tiles are given as 0-360
    df["two_scale"] = df["large_wl"] != df["small_wl"]
    df = df.reset_index(drop=True)
    df.insert(0, "tile", df.index)
    return df


if __name__ == "__main__":
    t = load_tiles()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT, index=False)
    print(t.shape, "tiles;", t["two_scale"].sum(), "with two dune scales")
    print(t.groupby("field").agg(n=("tile", "size"), two=("two_scale", "sum")).to_string())
    print(t.loc[t.two_scale, ["large_wl", "small_wl"]].describe().round(0))
