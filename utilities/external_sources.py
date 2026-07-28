"""
Non-ISMIP6 published ice-sheet simulations, downloaded directly from their
original public archives and processed here into the same shape as
utilities/data_loader.py's ismip6 dataframes: Year, Cumulative ice sheet mass
change (Gt), Group, Model, Exp, IS (+ a few extra classification columns).
This lets them be toggled on/off alongside ISMIP6 in the interactive
rate-comparison figure (see plot_interactive_rate_comparison's
`extra_sources` parameter).

All three sources were selected via a literature search (The Cryosphere/
Nature/GRL) for papers that (a) simulate future ice-sheet change, (b) also
include a historical/near-term run overlapping (or close to) the IMBIE
period, and (c) have a publicly archived dataset of annual-resolution mass
values -- not just plots in the paper.

Every raw file this module reads is downloaded on demand from its original
archive (NIRD Research Data Archive / Zenodo / Arctic Data Center) and
cached locally in CACHE_DIR, so a second run doesn't re-fetch it -- delete
CACHE_DIR (or the specific file inside it) to force a fresh download.
Nothing here reads from a pre-built/derived CSV: every call to
load_rahlves2025_gis()/load_coulon2024_ais()/load_aschwanden2022_gis()
re-derives the returned dataframe from the cached raw archive files, so the
processing steps below are always what actually produced the data you're
looking at, not a stale snapshot.
"""

import os
import re
import urllib.request

import h5py
import numpy as np
import pandas as pd
import scipy.io
import xarray as xr

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "external_sources_cache")


def _download(url, cache_name, min_expected_bytes=1):
    """Downloads `url` to CACHE_DIR/cache_name if not already cached, and
    returns the local path. This is the only function in this module that
    talks to the network -- everything else operates on files it returns.

    No checksum verification against the archive's own hashes (both
    archives publish per-file MD5s in their metadata; this only checks the
    download didn't obviously truncate/fail). If a download is interrupted,
    delete the partial file from CACHE_DIR and re-run.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    local_path = os.path.join(CACHE_DIR, cache_name)
    if os.path.exists(local_path) and os.path.getsize(local_path) >= min_expected_bytes:
        return local_path
    print(f"external_sources: downloading {cache_name} ...")
    urllib.request.urlretrieve(url, local_path)
    return local_path


# ═════════════════════════════════════════════════════════════════════════
# Rahlves et al. (2025), "Historically consistent mass loss projections of
# the Greenland ice sheet", The Cryosphere 19, 1205,
# https://doi.org/10.5194/tc-19-1205-2025.
#
# CISM (Weertman power-law sliding, inverse-calibrated basal friction),
# ERA5-forced historical run (tracks IMBIE within observed uncertainty) and
# ESM-forced runs, projected under SSP1-2.6/2-4.5/5-8.5 (RCP8.5 for the
# CMIP5-forced member) to 2100/2101, at three ocean-sensitivity tiers
# (Low/Medium/High, mapped onto the same `ocean_sensitivity` field ISMIP6's
# own GIS experiments use).
#
# Archive: NIRD Research Data Archive, https://doi.org/10.11582/2024.00128
# (767 GB total -- per-run `scalars.nc` files, the same domain-integrated-
# scalar convention ISMIP6 itself uses, following individual per-run
# directories). We don't download the whole archive: Step 1 below fetches
# just its own table-of-contents CSV (a file the archive publishes listing
# every file's name/URL/size) and filters it in code to the small subset
# actually used, so the file selection is visible here rather than baked
# into a hardcoded URL list.
# ═════════════════════════════════════════════════════════════════════════

RAHLVES2025_TOC_URL = (
    "https://data.archive.sigma2.no/dataset/ced83cf0-cf56-4289-8e95-b29a32395516/"
    "download/table_of_contents_10.11582_2024.00128.csv"
)

R_TIER_TO_OCEAN_SENSITIVITY = {"Low": "Low", "Med": "Medium", "High": "High"}
ESM_SCENARIO_MAP = {
    "hist": "Historical", "proj_ctrl": "Control", "proj_rcp85": "RCP 8.5",
    "proj_ssp126": "SSP1-2.6", "proj_ssp245": "SSP2-4.5", "proj_ssp585": "SSP5-8.5",
}
ERA5_F_CODE_TO_SCENARIO = {"26": "SSP1-2.6", "45": "SSP2-4.5", "85": "SSP5-8.5"}


def _rahlves2025_select_scalars_urls():
    """
    Step 1 of loading Rahlves et al. (2025): download the archive's own
    table of contents (every file it contains, with a download URL each --
    ~2.5 MB, one row per file across 767 GB of gridded model output) and
    filter it down to just the per-run scalar (domain-integrated) time
    series actually needed here.

    Filter applied, and why:
    - Only `scalars.nc` files (not the gridded 2D fields alongside them).
    - Only the "hist" (historical, ends 2015) and forward-projection runs,
      not the spin-up/relaxation runs that precede "hist" in each pipeline.
    - Only ONE resolution per initialization branch: 04 km for the
      ERA5-initialized branch (which also has 08/16 km resolution-
      sensitivity variants of the *same* runs -- including those too would
      triple-count identical parameter combinations at different grid
      resolutions -- 04 km is the archive's native/finest resolution for
      this branch) and 04 km for the ESM-initialized branch (its only
      resolution in this archive).
    - Excludes three one-off "proj_ctrl_o25"/"_o50"/"_o75" runs (one per
      ocean-sensitivity tier, each at a different ocean-forcing offset) that
      exist ONLY at 04 km, not the paper's regular Control scenario (that's
      the separate, symmetric "..._ctrl_proj" run present at all three
      tiers and both resolutions, and IS kept) -- these three don't fit the
      filename convention _rahlves2025_parse_filename decodes for every
      other run, and being an inconsistent one-per-tier sensitivity check
      rather than a systematic combination, they don't look like part of
      the paper's core analyzed ensemble.

    Returns {filename: download_url} for every file selected.
    """
    toc_path = _download(RAHLVES2025_TOC_URL, "rahlves2025_table_of_contents.csv", min_expected_bytes=10_000)
    toc = pd.read_csv(toc_path, sep="|", engine="python", skipinitialspace=True)
    toc.columns = [c.strip() for c in toc.columns]
    toc["http_url"] = toc["http_url"].str.strip()
    toc["filename"] = toc["filename"].str.strip()

    is_scalars = toc["filename"].str.endswith("scalars.nc")
    is_spinup_relax = toc["filename"].str.contains("spinup|relax", regex=True)
    is_ctrl_offset_outlier = toc["filename"].str.contains("proj_ctrl_o")
    is_era5_04km = toc["filename"].str.contains("ERA5_init/04km/")
    is_esm_04km = toc["filename"].str.contains("ESM_init/")  # only resolution present

    selected = toc[is_scalars & ~is_spinup_relax & ~is_ctrl_offset_outlier & (is_era5_04km | is_esm_04km)]
    # Every file is literally named "scalars.nc" -- the run identity lives
    # entirely in the directory path, so key on the full relative filename
    # (not os.path.basename, which would collapse all 162 to one key).
    return dict(zip(selected["filename"], selected["http_url"]))


def _rahlves2025_download_scalars():
    """Step 2: download (if not already cached) every scalars.nc file
    _rahlves2025_select_scalars_urls() selected -- 162 files, ~18-21 KB
    each, ~3.6 MB total. Returns {filename: local_path}."""
    urls_by_name = _rahlves2025_select_scalars_urls()
    local_paths = {}
    for url in urls_by_name.values():
        # Path components after the dataset root double as a unique, safe
        # local filename (same run can share a bare "scalars.nc" name).
        rel = url.split("Greenland_ice_sheet_projections/")[-1]
        cache_name = "rahlves2025_" + rel.replace("/", "__")
        local_paths[cache_name] = _download(url, cache_name, min_expected_bytes=1000)
    return local_paths


def _rahlves2025_parse_filename(cache_name):
    """Decodes one scalars.nc file's run identity from its path, following
    the archive's own directory-naming convention (see the module-level
    comment above for what each part means)."""
    core = cache_name.replace("rahlves2025_Greenland_ensemble_", "").replace("__scalars.nc", "")
    parts = core.split("__")
    init_group = parts[0]  # "ERA5_init" or "ESM_init"
    if init_group == "ESM_init":
        gcm = parts[1].split("_MAR312_", 1)[1]
        r_tier = re.search(r"_R(low|med|high)", parts[2], re.I).group(1).capitalize()
        runtype = parts[3]
        return {
            "init_group": "ESM", "gcm": gcm, "r_tier": r_tier, "runtype": runtype,
            "scenario": ESM_SCENARIO_MAP[runtype], "is_hist": runtype == "hist",
        }
    else:
        r_tier = re.search(r"_R(low|med|high)", parts[3], re.I).group(1).capitalize()
        runtype = parts[4]
        if runtype == "hist":
            return {"init_group": "ERA5", "gcm": "ERA5", "r_tier": r_tier, "runtype": runtype,
                    "scenario": "Historical", "is_hist": True}
        elif runtype.endswith("ctrl_proj"):
            return {"init_group": "ERA5", "gcm": "ERA5", "r_tier": r_tier, "runtype": runtype,
                    "scenario": "Control", "is_hist": False}
        else:
            m = re.search(r"_m(\d+)_r(\d+)_f(\d+)_o(\d+)", runtype)
            m_code, _, f_code, _ = m.groups()
            return {"init_group": "ERA5", "gcm": f"GCM-{m_code} (MAR-downscaled)", "r_tier": r_tier,
                    "runtype": runtype, "scenario": ERA5_F_CODE_TO_SCENARIO[f_code], "is_hist": False}


def load_rahlves2025_gis():
    """
    Loads Rahlves et al. (2025)'s Greenland/CISM ensemble as a dataframe
    shaped like ismip6_gis: Year, Cumulative ice sheet mass change (Gt),
    Group, Model, Exp, IS -- plus climate_model/scenario/protocol/
    ocean_sensitivity columns for classification (see
    utilities/external_sources.exp_meta_from_df).

    Step 1-2 (above): download the archive's table of contents, select and
    download the 162 relevant scalars.nc files.
    Step 3: pair each projection run with its matching historical (1961-
    2015) lead-in from the same init-group/R-tier (ESM_init) or init-group/
    GCM/R-tier (ERA5_init), and concatenate them into one continuous
    trajectory -- ISMIP6's own experiments are structured the same way
    (historical + projection concatenated), so this keeps the two directly
    comparable.
    Step 4: read each file's `imass_above_flotation` (ice mass above
    flotation, kg -- the CISM scalar variable ISMIP6 itself reports)
    and convert to Gt.
    """
    local_paths = _rahlves2025_download_scalars()
    meta_by_file = {name: _rahlves2025_parse_filename(name) for name in local_paths}

    hist_by_key = {}
    for name, m in meta_by_file.items():
        if m["is_hist"]:
            key = (m["init_group"], m["r_tier"]) if m["init_group"] == "ERA5" else (m["init_group"], m["gcm"], m["r_tier"])
            hist_by_key[key] = name

    rows = []
    exp_meta_rows = []
    for name, m in meta_by_file.items():
        if m["is_hist"]:
            continue  # merged into each projection below, not its own standalone Exp
        key = (m["init_group"], m["r_tier"]) if m["init_group"] == "ERA5" else (m["init_group"], m["gcm"], m["r_tier"])
        hist_name = hist_by_key.get(key)

        years, mass_gt = [], []
        if hist_name is not None:
            ds_h = xr.open_dataset(local_paths[hist_name])
            years.extend([t.year for t in ds_h.time.values])
            mass_gt.extend((ds_h["imass_above_flotation"].values / 1e12).tolist())
            ds_h.close()

        ds_p = xr.open_dataset(local_paths[name])
        p_years = [t.year for t in ds_p.time.values]
        p_mass = (ds_p["imass_above_flotation"].values / 1e12).tolist()
        ds_p.close()
        if years and p_years and p_years[0] == years[-1]:
            p_years, p_mass = p_years[1:], p_mass[1:]  # avoid a duplicate year at the hist/proj seam
        years.extend(p_years)
        mass_gt.extend(p_mass)

        exp_code = f"{m['init_group']}_{m['r_tier']}_{m['gcm']}_{m['scenario']}".replace(" ", "").replace("/", "-")
        rows.append(pd.DataFrame({
            "Year": years, "Cumulative ice sheet mass change (Gt)": mass_gt,
            "Group": "Rahlves2025", "Model": "CISM", "Exp": exp_code, "IS": "GIS",
        }))
        exp_meta_rows.append({
            "Exp": exp_code, "climate_model": m["gcm"], "scenario": m["scenario"],
            "protocol": "Not applicable", "ocean_sensitivity": R_TIER_TO_OCEAN_SENSITIVITY[m["r_tier"]],
        })

    df = pd.concat(rows, ignore_index=True)
    exp_meta_df = pd.DataFrame(exp_meta_rows)
    return df.merge(exp_meta_df, on="Exp", how="left")


# ═════════════════════════════════════════════════════════════════════════
# Coulon, Klose, Kittel, Edwards, Turner, Winkelmann & Pattyn (2024),
# "Disentangling the drivers of future Antarctic ice loss with a
# historically calibrated ice-sheet model", The Cryosphere 18, 653,
# https://doi.org/10.5194/tc-18-653-2024.
#
# Kori-ULB (f.ETISh), 100-member Latin-hypercube parameter ensemble,
# historical run 1950-2014 forced by NorESM1-M (CMIP5), Bayesian-calibrated
# against IMBIE regional mass balance (Otosaka et al., 2023), extended
# through SSP585_KEEP_THROUGH_YEAR using each member's own SSP5-8.5
# continuation (2015-3014 in the archive) from the same record. This uses
# the full unweighted 100-member ensemble (every parameter combination
# sampled, not just the Bayesian-calibrated best estimate) -- every mass
# trajectory the archive provides for the historical period, nothing
# pre-selected.
#
# Archive: Coulon et al. (2023), Zenodo, https://doi.org/10.5281/zenodo.8398771
# (confirmed directly from the paper's own "Code and data availability"
# section -- a different, unrelated Kori-ULB Zenodo record surfaced first
# when searching by model name alone, so don't trust a DOI found any other
# way for this dataset).
# ═════════════════════════════════════════════════════════════════════════

# Zenodo's concept DOI (8398771, covers all versions) resolves to a specific
# version record via the API; this is that resolved version's file API base.
COULON2024_FILES_BASE = "https://zenodo.org/api/records/10812218/files"

MELT_PARAM_MAP = {
    1: "Quadratic-local (Antarctic slope)", 2: "PICO", 3: "Plume",
    4: "ISMIP6 Nonlocal quadratic", 5: "ISMIP6 Nonlocal quadratic (slope-dependent)",
}
CLIMATE_MODEL_MAP = {1: "MRI-ESM2-0", 2: "UKESM1-0-LL", 3: "CESM2-WACCM", 4: "IPSL-CM6A-LR"}
GT_PER_M_SLE = 362.5  # Gt per meter of sea-level-equivalent
SSP585_KEEP_THROUGH_YEAR = 2030  # trim the 2015-3014 continuation here; see _coulon2024_download_ssp585


def _coulon2024_download_files():
    """
    Step 1: download the two archive files actually needed (out of ~46 in
    the full record, most of which are gridded 2D fields not used here):
    - HIST_ENSEMBLE_DATA.mat (~470 MB): historical (1950-2014) ensemble,
      including SLC_ensemble (cumulative sea-level contribution per member).
    - LHSensemble.mat (~11 KB): the 100-member Latin-hypercube parameter
      table (documented in the archive's own README.txt), giving each
      member's basal-melt parameterization and CMIP6 forcing GCM.
    Both are MATLAB files -- HIST_ENSEMBLE_DATA.mat is v7.3 (HDF5-backed,
    needs h5py), LHSensemble.mat is an older format (needs scipy.io).
    """
    hist_path = _download(
        f"{COULON2024_FILES_BASE}/HIST_ENSEMBLE_DATA.mat/content",
        "coulon2024_HIST_ENSEMBLE_DATA.mat", min_expected_bytes=100_000_000,
    )
    lhs_path = _download(
        f"{COULON2024_FILES_BASE}/LHSensemble.mat/content",
        "coulon2024_LHSensemble.mat", min_expected_bytes=1000,
    )
    return hist_path, lhs_path


def _coulon2024_download_ssp585():
    """
    Step 4: download SSP585_ENSEMBLE_DATA.mat, the SSP5-8.5 continuation
    (2015-3014) of the same 100-member ensemble, per the archive's own
    README.txt -- used here to extend each member's historical (1950-2014)
    trajectory a bit further into the projection period.

    A partial/byte-range download (fetching only the leading portion of the
    file, since its 1.2 GB is almost entirely the gridded H_ensemble/
    MASK_ensemble/Runoff_ensemble fields this module doesn't use, not the
    tiny SLC_ensemble time series it actually needs) was checked first:
        curl -sI -H "Range: bytes=0-1023" \\
            https://zenodo.org/api/records/10812218/files/SSP585_ENSEMBLE_DATA.mat/content
    returns a plain "200 OK" with the full Content-Length (no
    "Accept-Ranges"/"Content-Range" headers, i.e. a 206 Partial Content
    never comes back) -- Zenodo's file-serving endpoint doesn't support
    HTTP Range requests, so there is no way to fetch only part of this file
    over HTTP. Falling back to a full download, as a result. `time`/
    `SLC_ensemble` are read with h5py in Step 5 below, which only pulls
    those two (small) datasets' bytes out of the local file, not the whole
    thing into memory; anything past SSP585_KEEP_THROUGH_YEAR is discarded
    immediately after that read.
    """
    return _download(
        f"{COULON2024_FILES_BASE}/SSP585_ENSEMBLE_DATA.mat/content",
        "coulon2024_SSP585_ENSEMBLE_DATA.mat", min_expected_bytes=1_000_000_000,
    )


def load_coulon2024_ais():
    """
    Loads Coulon et al. (2024)'s Antarctic/Kori-ULB 100-member ensemble as
    a dataframe shaped like ismip6_ais: Year, Cumulative ice sheet mass
    change (Gt), Group, Model, Exp, IS -- plus climate_model/scenario/
    protocol/basal_melt_param columns for classification.

    Step 1 (above): download HIST_ENSEMBLE_DATA.mat and LHSensemble.mat.
    Step 2: read HIST_ENSEMBLE_DATA.mat's SLC_ensemble (cumulative
    sea-level contribution, m, per member per year, 1950-2014).
    Step 3: read LHSensemble.mat's basal-melt-parameterization (column 7)
    and CMIP6-forcing-GCM (column 9) columns, decoded per the archive's own
    README.txt, to classify each of the 100 members.
    Step 4 (see _coulon2024_download_ssp585): download
    SSP585_ENSEMBLE_DATA.mat and read its own SLC_ensemble too, trimmed to
    SSP585_KEEP_THROUGH_YEAR. Its SLC_ensemble is NOT a continuation of
    HIST_ENSEMBLE_DATA's cumulative integral -- inspecting the raw values
    confirms it resets to 0 at 2015 (it's its own standalone 2015-3014
    simulation, re-zeroed at its own start year, same as HIST_ENSEMBLE_DATA
    is re-zeroed at 1950) -- so each member's SSP5-8.5 segment is shifted by
    that member's own HIST_ENSEMBLE_DATA endpoint (its cumulative SLC as of
    2014) before concatenating, restoring one continuous integral. Skipping
    this shift was tried first and produces an unphysical jump at the
    2014/2015 seam (e.g. member000 goes from -13914 Gt in 2014 to -0 Gt in
    2015), which shows up as wildly inflated 2010-2020 rates (mean ~980
    Gt/yr, max ~15000 Gt/yr) -- clearly wrong for a 100-member ensemble
    whose historical rates individually sit in the tens-to-low-hundreds of
    Gt/yr.
    Step 5: concatenate each member's historical and shifted SSP5-8.5
    SLC_ensemble into one continuous 1950-SSP585_KEEP_THROUGH_YEAR
    trajectory and convert to Gt (1 m SLE = 362.5 Gt; positive SLC = mass
    loss, hence the sign flip) -- the same historical+projection
    concatenation load_rahlves2025_gis() does, so the two stay directly
    comparable. Every member's `scenario` is labelled "SSP5-8.5" (the
    forcing every returned trajectory now runs on after 2014), matching how
    load_rahlves2025_gis() labels its own merged runs by their forward
    scenario rather than "Historical".
    """
    hist_path, lhs_path = _coulon2024_download_files()
    ssp585_path = _coulon2024_download_ssp585()

    with h5py.File(hist_path, "r") as f:
        hist_time = np.array(f["time"]).squeeze()  # (65,) years 1950-2014
        hist_slc = np.array(f["SLC_ensemble"])  # (65, 100) meters SLE, cumulative from 1950

    with h5py.File(ssp585_path, "r") as f:
        ssp_time = np.array(f["time"]).squeeze()  # years 2015-3014
        ssp_slc = np.array(f["SLC_ensemble"])  # (n_years, 100) meters SLE, re-zeroed at 2015
    keep = ssp_time <= SSP585_KEEP_THROUGH_YEAR
    ssp_time, ssp_slc = ssp_time[keep], ssp_slc[keep]

    assert ssp_time[0] == hist_time[-1] + 1, (
        f"expected the SSP5-8.5 file to pick up the year after HIST_ENSEMBLE_DATA "
        f"ends, but HIST ends {hist_time[-1]} and SSP585 starts {ssp_time[0]}"
    )
    ssp_slc = ssp_slc + hist_slc[-1, :]  # re-anchor each member onto its own HIST endpoint

    lhs_mat = scipy.io.loadmat(lhs_path)
    lhval = np.array(lhs_mat["LHval"])
    if lhval.shape[0] == 9:
        lhval = lhval.T  # README documents this as 100x9; transpose if MATLAB stored it column-major
    assert lhval.shape == (100, 9), f"unexpected LHSensemble.mat LHval shape {lhval.shape}"
    melt_param_col = lhval[:, 6]  # 7th column (1-indexed)
    climate_model_col = lhval[:, 8]  # 9th column

    time = np.concatenate([hist_time, ssp_time])
    rows = []
    exp_meta_rows = []
    for member in range(100):
        slc = np.concatenate([hist_slc[:, member], ssp_slc[:, member]])
        cum_gt = -slc * 1000 * GT_PER_M_SLE  # m SLE -> mm -> Gt, sign-flipped
        exp_code = f"member{member:03d}"
        rows.append(pd.DataFrame({
            "Year": time.astype(int), "Cumulative ice sheet mass change (Gt)": cum_gt,
            "Group": "Coulon2024", "Model": "Kori-ULB", "Exp": exp_code, "IS": "AIS",
        }))
        exp_meta_rows.append({
            "Exp": exp_code,
            "climate_model": CLIMATE_MODEL_MAP.get(int(round(climate_model_col[member])), "Unknown"),
            "scenario": "SSP5-8.5", "protocol": "Not applicable",
            "basal_melt_param": MELT_PARAM_MAP.get(int(round(melt_param_col[member])), "Unknown"),
        })

    df = pd.concat(rows, ignore_index=True)
    exp_meta_df = pd.DataFrame(exp_meta_rows)
    return df.merge(exp_meta_df, on="Exp", how="left")


# ═════════════════════════════════════════════════════════════════════════
# Aschwanden & Brinkerhoff (2022), "Calibrated Mass Loss Predictions for the
# Greenland Ice Sheet", Geophysical Research Letters 49, e2022GL099058,
# https://doi.org/10.1029/2022GL099058.
#
# PISM, forced by a temperature-index SMB model driven by RCP-scenario
# warming (not per-GCM downscaling the way Rahlves2025/Coulon2024 are, so
# there's no climate_model classification here beyond the RCP itself). The
# paper trains a PyTorch neural-network emulator on ~1000 true PISM runs,
# then draws large ensembles from it both before ("les", the prior --
# what's loaded here, per an explicit user request to use the uncalibrated
# ensemble rather than the Bayesian-calibrated posterior) and after ("mc")
# jointly conditioning on observed surface speed and observed cumulative
# mass loss.
#
# Archive: Aschwanden & Brinkerhoff (2022), Arctic Data Center,
# https://doi.org/10.18739/A2KW57K4R -- a plain Apache directory listing
# (https://arcticdata.io/data/10.18739/A2KW57K4R/pism_scalars/), no
# API/DOI-resolution archaeology needed the way Coulon2024's Zenodo record
# or Rahlves2025's NIRD table-of-contents did.
# ═════════════════════════════════════════════════════════════════════════

ASCHWANDEN2022_LES_URL = (
    "https://arcticdata.io/data/10.18739/A2KW57K4R/pism_scalars/"
    "aschwanden_et_al_2019_les_2008_norm.csv.gz"
)
ASCHWANDEN2022_RCP_TO_SCENARIO = {26: "RCP2.6", 45: "RCP4.5", 85: "RCP8.5"}
ASCHWANDEN2022_KEEP_THROUGH_YEAR = 2100  # trim the 2008-3007 file to the paper's own stated analysis horizon


def _aschwanden2022_download_les():
    """
    Step 1: download the archive's uncalibrated prior ensemble file
    (~56 MB gzipped CSV) -- one row per (Experiment, RCP, Year), zero-
    referenced at 2008 (its own "_norm" filename suffix).

    This is a plain CSV, unlike Coulon2024/Rahlves2025's MATLAB/NetCDF
    files, so no h5py/scipy.io/xarray parsing is needed below -- just
    pandas.
    """
    return _download(ASCHWANDEN2022_LES_URL, "aschwanden2022_les_2008_norm.csv.gz", min_expected_bytes=10_000_000)


def load_aschwanden2022_gis():
    """
    Loads Aschwanden & Brinkerhoff (2022)'s Greenland/PISM uncalibrated
    prior ensemble as a dataframe shaped like ismip6_gis: Year, Cumulative
    ice sheet mass change (Gt), Group, Model, Exp, IS -- plus
    climate_model/scenario/protocol columns for classification.

    Step 1 (above): download aschwanden_et_al_2019_les_2008_norm.csv.gz.
    Step 2: read the `Mass (Gt)` column directly (already the exact
    cumulative-mass-change quantity this module's other loaders have to
    derive via unit conversion) for each (Experiment, RCP) realization,
    trimmed to ASCHWANDEN2022_KEEP_THROUGH_YEAR (the file runs to year
    3007, like Coulon2024's SSP585 continuation, far past anything this
    figure's year-range slider needs).
    Step 3: label each realization's `scenario` from its RCP and mark
    `protocol` as "Uncalibrated (prior ensemble)" -- distinguishing it from
    the archive's separate Bayesian-calibrated "mc" (posterior) ensemble,
    which this loader does NOT use.

    Zero-referenced at 2008 rather than at a common absolute ice-sheet mass
    -- like every other source in this module, only rate differences over a
    selected year window are ever computed from this, so the absolute
    reference year doesn't matter, but it does mean a rate window starting
    before 2008 will have no Aschwanden2022 points.
    """
    les_path = _aschwanden2022_download_les()
    df = pd.read_csv(les_path)
    df = df[df["Year"] <= ASCHWANDEN2022_KEEP_THROUGH_YEAR].copy()

    rcp = df["RCP"].astype(int)
    exp = df["Experiment"].astype(int)
    df["Exp"] = [f"les{e:03d}_rcp{r}" for e, r in zip(exp, rcp)]
    df["Year"] = df["Year"].astype(int)
    df["Cumulative ice sheet mass change (Gt)"] = df["Mass (Gt)"]
    df["Group"] = "Aschwanden2022"
    df["Model"] = "PISM"
    df["IS"] = "GIS"
    df["climate_model"] = "Not applicable"
    df["scenario"] = rcp.map(ASCHWANDEN2022_RCP_TO_SCENARIO)
    df["protocol"] = "Uncalibrated (prior ensemble)"

    return df[[
        "Year", "Cumulative ice sheet mass change (Gt)", "Group", "Model", "Exp", "IS",
        "climate_model", "scenario", "protocol",
    ]]


def exp_meta_from_df(df, extra_cols):
    """Builds a get_exp_meta()-style {Exp: {...}} dict from one of this
    module's loaded dataframes, whose exp-level classification (scenario,
    climate_model, protocol, and either ocean_sensitivity or
    basal_melt_param) is already attached as columns -- unlike ISMIP6's
    hand-transcribed gis_exp_meta/ais_exp_meta, there's no separate lookup
    table to maintain here."""
    cols = ["Exp", "climate_model", "scenario", "protocol"] + extra_cols
    unique = df[cols].drop_duplicates("Exp").set_index("Exp")
    return unique.to_dict(orient="index")
