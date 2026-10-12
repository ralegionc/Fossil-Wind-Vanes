"""Find PMIP monthly surface-wind files (uas, vas) on ESGF and write a manifest.

For each model and experiment we keep one ensemble member and one dataset
version, and only the files that cover the last 50 simulated years. The
manifest lists HTTP and OPeNDAP URLs so the next step can subset remotely.
"""
from collections import defaultdict
import json
from pathlib import Path
import re

import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed" / "pmip_manifest.json"
SEARCH = "https://esgf.ceda.ac.uk/esg-search/search"
YEARS = 50

CMIP6 = ["AWI-ESM-1-1-LR", "INM-CM4-8", "MIROC-ES2L", "MPI-ESM1-2-LR"]
CMIP5 = ["CNRM-CM5", "GISS-E2-R", "IPSL-CM5A-LR", "MIROC-ESM", "MPI-ESM-P", "MRI-CGCM3"]
EXPS = ["lgm", "piControl", "midHolocene"]
VARS = ["uas", "vas"]


def search_files(project, model, exp, var):
    if project == "CMIP6":
        q = dict(project="CMIP6", source_id=model, experiment_id=exp, variable_id=var, table_id="Amon")
    else:
        q = dict(project="CMIP5", model=model, experiment=exp, variable=var, time_frequency="mon",
                 cmor_table="Amon")
    q.update(type="File", distrib="true", limit=10000, format="application/solr+json",
             fields="id,title,url,size,data_node,dataset_id,version,member_id,ensemble")
    r = requests.get(SEARCH, params=q, timeout=120)
    r.raise_for_status()
    return r.json()["response"]["docs"]


def member_of(doc):
    m = doc.get("member_id") or doc.get("ensemble") or ["?"]
    return m[0] if isinstance(m, list) else m


def span(title):
    m = re.search(r"_(\d{4,6})(\d{2})?-(\d{4,6})(\d{2})?\.nc$", title)
    if not m:
        return None
    a, b = m.group(1), m.group(3)
    return int(a[:-2] if len(a) == 6 else a), int(b[:-2] if len(b) == 6 else b)


# GISS-E2-R ran its LGM with physics version p150, built on the p141 pre-industrial
# configuration, so its LGM change is taken against the p141 control.
EXTRA = {("GISS-E2-R", "piControl"): ("piControl_p141", "r1i1p141")}


def pick(docs, prefer=None):
    """One member (prefer r1i1p1[f1]), newest dataset version, files for the last YEARS years."""
    if not docs:
        return None
    by_member = defaultdict(list)
    for d in docs:
        by_member[member_of(d)].append(d)
    order = ([prefer] if prefer else []) + ["r1i1p1f1", "r1i1p1"]
    member = next((m for m in order if m in by_member), sorted(by_member)[0])
    docs = by_member[member]
    # dataset_id looks like <dataset>.v<version>|<node>; keep the newest version
    def ver(d):
        return d["dataset_id"].split("|")[0].rsplit(".", 1)[-1]
    newest = max(ver(d) for d in docs)
    docs = [d for d in docs if ver(d) == newest]
    files = {}
    for d in docs:  # same file can sit on several nodes; collect all URLs
        f = files.setdefault(d["title"], {"title": d["title"], "size": d.get("size", 0), "http": [], "opendap": []})
        for u in d["url"]:
            url, _, kind = u.split("|")
            if kind == "HTTPServer":
                f["http"].append(url)
            elif kind == "OPENDAP":
                f["opendap"].append(url.removesuffix(".html"))
    spans = {t: span(t) for t in files}
    end = max(s[1] for s in spans.values() if s)
    keep = [files[t] for t, s in spans.items() if s and s[1] > end - YEARS]
    return {"member": member, "version": newest, "end_year": end,
            "files": sorted(keep, key=lambda f: f["title"])}


def main():
    manifest = {}
    for project, models in (("CMIP6", CMIP6), ("CMIP5", CMIP5)):
        for model in models:
            for exp in EXPS:
                for var in VARS:
                    docs = search_files(project, model, exp, var)
                    jobs = [(f"{model}|{exp}|{var}", None)]
                    if (model, exp) in EXTRA:
                        name, member = EXTRA[(model, exp)]
                        jobs.append((f"{model}|{name}|{var}", member))
                    for key, prefer in jobs:
                        sel = pick(docs, prefer)
                        manifest[key] = sel
                        if sel:
                            mb = sum(f["size"] for f in sel["files"]) / 1e6
                            print(f"{key:40s} {sel['member']:10s} {len(sel['files']):3d} files {mb:7.1f} MB end {sel['end_year']}")
                        else:
                            print(f"{key:40s} not found")
    OUT.write_text(json.dumps(manifest, indent=1))


if __name__ == "__main__":
    main()
