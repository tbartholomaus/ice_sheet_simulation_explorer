# Metadata review 3: Goelzer et al. (2025), PROTECT Greenland ensemble

**Paper:** Goelzer, H., et al. (2025). Extending the range and reach of physically-based Greenland ice sheet sea-level projections. *The Cryosphere* 19, 6887–6906. https://doi.org/10.5194/tc-19-6887-2025 (CC BY 4.0).

**Data:** `dash_app/data/external_sources_goelzer2025_gis.csv.gz`: 755 runs, built by `utilities/external_sources.py` `load_goelzer2025_gis()` from the NIRD archive (https://doi.org/10.11582/2025.lf9m2wd0). The build keeps only the 2015–2100 scenarios and, for NORCE, only the finest-resolution copy of each experiment.

**Rules applied:** those of reviews 1 and 2 (`encoding_rules.csv`).

---

## 1. Summary

| # | Severity | Finding |
|---|---|---|
| 1 | **High** | **Each lab ran several distinct model configurations, but the data and metadata collapse each lab to one model** (Group "Goelzer2025", Model = Elmer/Ice, IMAU-ICE, CISM or GISM). Table B1 lists the configurations and what distinguishes them. Several differ in exactly the characteristics the app analyses: ElmerIce2 vs ElmerIce3 (linear vs Weertman sliding); IMAUICE1–8 (10–40 km resolution, three spin-up forcings, Zoet–Iverson vs regularized Coulomb sliding); NORCE CISM (2–16 km, MARv3.9 vs MARv3.12 spin-up, historical retreat sensitivity); GISMHOMv1 vs GISMSIAv1 (higher-order vs shallow-ice). This conflicts with the configuration-independence rule. Decision G1. |
| 2 | **High** | The current (unreviewed) encodings are wrong for most configurations under the agreed classes; §3 gives the reviewed values. |
| 3 | Medium | *(Resolved by G2.)* **NORCE runs were reduced to the finest resolution** per experiment (a 2026-08-03 choice), so NORCE's 8 and 16 km reruns are missing, while all of IMAU's resolution variants (10–40 km) are kept. Decision G2. |
| 4 | Medium | **Retreat sensitivity has five levels** (p05, p25, p50, p75, p95 of the retreat parameter κ), where ISMIP6-GIS has three (25th/50th/75th = Low/Medium/High). Decision G3. |
| 5 | Medium | **The regional climate model is a forcing choice** here (MARv3.9, MARv3.12, RACMO2.3p2, HIRHAM5; Table 1) and isn't encoded anywhere. Decision G4. |
| 6 | Low | GCM names: "ACCESS1.3" and "CSIRO-Mk3.6.0" map to the CMIP identifiers ACCESS1-3 and CSIRO-Mk3-6-0 (rule D5). "CESM2-Leo" and "UKESM1-0-LL-Robin" are pre-CMIP6 ensemble members (Table 1 footnote a), so they are distinct from CESM2 / UKESM1-0-LL. Decision G5. |
| 7 | Low | Unlike ISMIP6, the archive values are not control-corrected. That is consistent with the paper: historical retreat forcing makes the pre-2015 run "a real historical simulation", so ISMIP6's control-run removal "is therefore not needed here" (Sect. 2.5.1). No change. |

---

## 2. Experiments

PROTECT has no single experiment table. Each experiment is a combination of [Global model, GHG scenario, Regional model, Retreat sensitivity] (Sect. 2.5.4; Fig. 4). Table 1 lists which ESM × scenario × RCM forcings exist. The archive encodes each experiment as `<GCM>_<scenario>_<RCM>_<retreat percentile>`, which the loader parses.

- **Scenarios:** SSP1-2.6, SSP2-4.5 and SSP5-8.5 (CMIP6), plus RCP2.6 and RCP8.5 (CMIP5). The current encoding is correct.
- **Ocean forcing:** the ISMIP6 standard retreat parameterization for every run ("the standard ISMIP6 protocol is applied", Sects. 2.4.2–2.4.3, 2.5.2) → **Standard**.
- **Retreat sensitivity:** p05–p95 of κ (Sect. 2.5.1; Fig. 6). See G3.
- **Regional model:** MARv3.9 (produced for ISMIP6), MARv3.12, RACMO2.3p2, HIRHAM5. See G4.

---

## 3. Model configurations (Table 2, Table B1, Sect. 2.4)

### 3a. Sources

| Key | Source | Used for |
|---|---|---|
| G25 | Goelzer et al. 2025, Table 2, Table B1, Sects. 2.4–2.5 | All configurations |
| B23 | Berends et al. 2023, TC 17, 1585 (cited by G25 for the IMAU inversion) | IMAU-ICE Zoet–Iverson sliding law; inversion against thickness and velocity |
| B22 | Berends et al. 2022, GMD 15, 5667 (IMAU-ICE v2.0, cited by G25) | IMAU-ICE regularized Coulomb law (τb = τc (|u|/u0)^q, q = 0.3, τc = tan φ × effective overburden) |
| F15 | Fürst et al. 2015, TC 9, 1039 (cited by G25) | GISM sliding: vb = SBL (AS/H) τb³ |
| G20 | Goelzer et al. 2020 App. A10 (review 2) | CISM: Schoof (2005)-based sliding, same wording as G25 Sect. 2.4.4 |

### 3b. Proposed encoding (one row per configuration, if G1 is accepted)

| Configuration (Table B1) | Stress balance | Res (km) → bin | Initialization → class | Sliding law → class |
|---|---|---|---|---|
| IGE ElmerIce2 | SSA → Shallow | 1–6 → ≤5 | Inverse control method on velocities (Gillet-Chaulet et al. 2012) + 20-yr relaxation → Data assimilation | "linear friction law" (G25 Sect. 2.4.3; Table B1 "Linear") → **Linear** |
| IGE ElmerIce3 | SSA → Shallow | 1–6 → ≤5 | as ElmerIce2 → Data assimilation | "Weertman friction law with m = 1/3" → **Power law (Weertman)** |
| IMAU IMAUICE1 / 2 / 3 / 5 | SIA-SSA → Shallow | 40 / 30 / 20 / 10 → >20 / >20 / >10–20 / >5–10 | Basal inversion of till friction angle against observed thickness and velocity (B23), then paleo spin-up from 120 ka (PMIP3 forcing) → **Hybrid (thickness-nudged)** | Zoet–Iverson, τb = N tan φ (|u|/(|u|+u0))^(1/p), N = ice overburden (B23 Eq. 1) → **Coulomb-limited** |
| IMAU IMAUICE6 / 7 | SIA-SSA → Shallow | 20 → >10–20 | as above, HadCM3 / CCSM spin-up forcing → Hybrid | Zoet–Iverson → Coulomb-limited |
| IMAU IMAUICE8 | SIA-SSA → Shallow | 20 → >10–20 | as above, PMIP3 → Hybrid | Regularized Coulomb (Table B1 "RC"; B22: τc = tan φ (ρi g H − pw)) → Coulomb-limited |
| NORCE CISM (each variant) | DIVA (depth-integrated higher-order, Goldberg 2011) → Higher-order | 2 / 4 / 8 / 16 → ≤5 / ≤5 / >5–10 / >10–20 | 5000-yr spin-up nudged toward present-day thickness by adjusting friction coefficients → Hybrid | "basal-sliding law based on Schoof (2005)"; "no dependence ... on basal temperature or water pressure" → Coulomb-limited |
| VUB GISMHOMv1 | HO → Higher-order | 5 → ≤5 | Assimilation of observed thickness (Le clec'h et al. 2019) → Hybrid | Weertman, vb ∝ (AS/H) τb³ (F15) → Budd |
| VUB GISMSIAv1 | SIA → Shallow | 5 → ≤5 | as GISMHOMv1 → Hybrid | as GISMHOMv1 → Budd |

NORCE configurations in the bundled data, decoded from the variant name and Table B1 (resolution; MAR version for spin-up and historical run; historical retreat percentile; "c" = retreat sensitivity consistent between historical and projection): CISM02-MAR39-p50; CISM04-MAR312-p25, -p50, -p75; CISM04-MAR39-p50; CISM04c-MAR39-p05/-p25/-p50/-p75/-p95; CISM08-MAR39-p50; CISM16c-MAR312-p05/-p25/-p50/-p75/-p95.

Notes:
- IMAU-ICE v2.0's default stress balance is DIVA (B22), but G25 Table 2 gives "SIA-SSA" for the PROTECT runs. The table is followed.
- IMAU-ICE Zoet–Iverson uses N = ice overburden, with no water pressure (B23).

---

## 4. Decisions needed

- ~~**G1. Split each lab into its configurations.**~~ **DECIDED:** yes. Group = "Goelzer2025-<lab>", Model = the Table B1 configuration; `ice_model` still names the code base (Elmer/Ice, IMAU-ICE, CISM, GISM), so grouping by ice sheet model shows four models.
- ~~**G2. NORCE resolution duplicates.**~~ **DECIDED:** restore them; every NORCE resolution is its own configuration.
- ~~**G3. Retreat sensitivity.**~~ **DECIDED:** p05/p25/p50/p75/p95 → Very low / Low / Medium / High / Very high (the two new values added to `encoding_rules.csv`).
- ~~**G4. Regional climate model.**~~ **DECIDED:** a new characteristic (grouping and ANOVA factor) across studies. ISMIP6-AIS: "None (AOGCM output used directly)" (Seroussi et al. 2020 Sect. 2.2); ISMIP6-GIS: MARv3.9 (Goelzer et al. 2020 Sect. 2); unreviewed studies: "Not yet reviewed".
- ~~**G5. Pre-CMIP6 members.**~~ **DECIDED:** CESM2-Leo and UKESM1-0-LL-Robin stay their own GCM values.


---

## 5. Implementation (2026-11-23)

- `utilities/external_sources.py` `load_goelzer2025_gis()`: Group = "Goelzer2025-<lab>", Model = the archive configuration, Exp = `<GCM>_<scenario>_<RCM>_<percentile>`, plus a `regional_climate_model` column. The NORCE finest-resolution selection is removed.
- `dash_app/data/external_sources_goelzer2025_gis.csv.gz` regenerated from the cached archive: **1270 runs, 39 configurations** (IGE 2, IMAU 7, NORCE 28, VUB 2). NORCE's per-configuration run counts match Table B1 for every 2015–2100 configuration. The 755 runs that were in the previous file are identical (values and forcing labels).
- `dash_app/metadata/Goelzer2025_models.csv` (39 reviewed rows) and `Goelzer2025_experiments.csv` (170 forcing combinations), with value maps for scenarios, GCM identifiers (ACCESS1.3 → ACCESS1-3, CSIRO-Mk3.6.0 → CSIRO-Mk3-6-0) and retreat percentiles.
- New characteristic `regional_climate_model` across studies (app grouping, ANOVA, hover; notebook). New ocean-sensitivity values Very low / Very high.
- Checks: 65 tests pass. Outside Goelzer 2025, the run table is unchanged in every cell (13,147 runs). The app's Greenland bias analysis runs with both Greenland ensembles (1435 runs, 11 factors), with no console errors.
