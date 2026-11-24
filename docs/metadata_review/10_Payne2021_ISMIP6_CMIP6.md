# Metadata review 10: Payne et al. (2021), ISMIP6 runs forced by CMIP6 (new study)

**Paper:** Payne, A. J., Nowicki, S., Abe-Ouchi, A., et al. (2021). Future sea level change under Coupled Model Intercomparison Project phase 5 and phase 6 scenarios from the Greenland and Antarctic ice sheets. *GRL* 48, e2020GL091741. https://doi.org/10.1029/2020GL091741 (open access, CC BY 4.0). The paper was read in full, including Tables 1–2. Citations to it below are "P21".

**Data:** https://doi.org/10.5281/zenodo.4498331 (CC BY 4.0), `CMIP5_CMIP6_Scalars_Paper.zip`:
- **Antarctica:** `computed_limnsw_minus_ctrl_proj_AIS_*_expB<n>.nc`, in Gt.
- **Greenland:** `scalars_mm_cr_GIS_*_expb<nn>.nc`. Per the ISMIP6-Greenland archive README, "suffix cr in filename for control removed".

Both are published with each configuration's control run subtracted.

**Decisions (user, 2026-11-24):**
1. **Add the control back,** so Payne matches the app's other, raw ISMIP6 runs.
2. **Map the pre-publication configuration names after verifying them.**
3. **Make it a separate study, weighted as ISMIP6.**

---

## 1. What was done and verified

| # | Item |
|---|---|
| 1 | **Experiments.** P21 Sect. 3: "(a) The CNRM-CM6-1 model run with scenario SSP585 ..., (b) CNRM-CM6-1 with SSP126 ..., and SSP585 with (c) UKESM1-0-LL, (d) CESM2, (e) CNRM-ESM2-1". Antarctic expB1–5 are the open framework and expB6–10 standard: Table 1 gives open and standard runs 1–5, and VUB, "standard 1–3" only, has exactly expB6–8. Greenland expb01–05 are standard and expb06–10 open (Table 2; BGC and UAF_PISM2, open only, have only b06–b10). |
| 2 | **Control added back.** See the note below the table. |
| 3 | **Renamed configurations.** The archive keeps four pre-publication names. The ISMIP6-Greenland analysis code's own name table (`github.com/ismip/ismip6-gris-analysis`, `define_set.m`, arrays `imod` → `papmod`) gives SICOPOLIS2 → SICOPOLIS1, SICOPOLIS3 → SICOPOLIS2, GRISLI2 → GRISLI and GISMHOMv1 → GISM. Each file's model constants (ρᵢ, ρ_w, time origin) also match its mapped control. Antarctic "SICOPOLIS1" is Seroussi 2020's only SICOPOLIS configuration (standard framework only, as in P21 Table 1). |
| 4 | **Counts.** Antarctica: 43 runs from 7 configurations (Table 1). Greenland: 69 runs from 16 configurations (Table 2; expb01, 14 runs = the 14 standard groups). |
| 5 | **Model metadata.** Every configuration is a Seroussi 2020 (AIS) or Goelzer 2020 (GIS) configuration, so each row copies that study's reviewed row (reviews 01, 02). A test checks that every class matches. |
| 6 | **Forcing.** The ESMs are CMIP6 identifiers already. Ocean sensitivity is **Medium**: P21 doesn't vary it, and the CMIP6 runs use ISMIP6's median settings. Regional climate model: Antarctica "None (AOGCM output used directly)"; Greenland **MARv3.9**. P21 says only "the MAR regional climate model of Greenland (Fettweis et al., 2013)"; the version is from Goelzer et al. 2025: "Forcing with MAR version 3.9 was produced for ISMIP6". |
| 7 | **App.** "Payne 2021 (ISMIP6 AIS, CMIP6)" and "Payne 2021 (ISMIP6 GIS, CMIP6)", in publication order after DeConto 2021, coloured #66a61e (CIEDE2000 ≥ 20.8 against every other study). Both are in `ISMIP6_STUDIES`, so each run weighs 1, and both are unchecked by default. The notebook adds them as extra sources. |

**How the control is added back (item 2).** The method is exactly that of the app's builders for Seroussi 2020 and Goelzer 2020 (`utilities/data_loader.py`, `remove_ctrl=False`):
- **Antarctica:** raw = P21 + `ctrl_proj_<open|std>` ivaf × ρᵢ, from the Seroussi 2020 archive. In that archive, raw − minus-control equals the control exactly (difference 0.0), and the P21 and control time axes are identical. The result is joined to `hist_<open|std>`.
- **Greenland:** raw = P21 + (control − control at 2015), from the Goelzer 2020 archive, joined to the historical run.
- **Both:** rebased to 0 at 2015.
- **Check:** the historical segments equal the app's existing runs for the same configurations to within 10⁻¹³ Gt.

**Caveat.** P21's Antarctic files give "limnsw" (land ice mass not displacing seawater, Gt), while the controls are ivaf × ρᵢ. These are the same quantity in ISMIP6, but the two archives compute them separately. Any difference would show as a small offset from 2016, which the 2015 rebasing doesn't remove. Proving the reconstruction exact would need P21's raw runs, which are on GHub behind a login.
