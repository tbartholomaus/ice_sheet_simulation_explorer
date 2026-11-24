# Metadata review 7: Rahlves et al. (2025), CISM Greenland ensemble

**Paper:** Rahlves, C., Goelzer, H., Born, A. and Langebroek, P. M. (2025). Historically consistent mass loss projections of the Greenland ice sheet. *The Cryosphere* 19, 1205–1220. https://doi.org/10.5194/tc-19-1205-2025 (open access). Sections 2–3, Table A1 and the figure captions were read. Citations to the paper below are "R25".

**Data:** NIRD Research Data Archive, https://doi.org/10.11582/2024.00128. The archive's table of contents and the 04 km ensemble READMEs were read.
- The loader takes the per-run `scalars.nc` files at 4 km for both branches.
- It joins each projection to its own historical run.
- In the app: 129 runs. ERA5-init has 17 per sensitivity tier (16 projections + 1 control); ESM-init has 26 per tier (16 projections + 10 controls).

**Rules applied:** those of reviews 1–6 (`encoding_rules.csv`).

**Also found during this review:** Goelzer 2025's retreat-percentile mapping was inverted (review 03, decision G3). It was corrected in the same pass; see the note under G3 in report 03.

---

## 1. Summary

| # | Severity | Finding |
|---|---|---|
| 1 | **High** | **The sliding law is mis-encoded** as "Weertman (unreviewed)". R25 Eq. 1: "τb = k^(−1/p) N^(q/p) ub^(1/p) ... N is the effective pressure ... p = 3 and q = 1". Friction depends on effective pressure, so cascade step 2 applies: **Budd**. |
| 2 | **High** | **The stress balance and initialization were never encoded.** CISM here "solv[es] a depth-integrated viscosity approximation of the Stokes equations (Goldberg, 2011)", which is DIVA, so **Higher-order**. The initialization is a 5000-year spin-up during which "basal friction parameters are calibrated to nudge the ice surface elevation towards observed conditions", then 1000 years of relaxation: **Hybrid (thickness-nudged)** (D1). |
| 3 | **High** | **The ERA5-init climate models are unnamed** ("GCM-06 (MAR-downscaled)" ... "GCM-15"). The archive README decodes them: `["ACCESS1.3"]="06" ["CNRM-CM6"]="07" ["UKESM1-CM6"]="08" ["CESM2"]="09" ["CNRM-ESM2"]="10" ["MPI-ESM1-2-HR"]="11" ["IPSL-CM6A-LR"]="12" ["NorESM2"]="13" ["CESM2-CMIP6"]="14" ["UKESM1-0-LL-CMIP6"]="15"`. The ten codes match Table A1's ten models in the same order (08 = UKESM1-0-LL-Robin, 09 = CESM2-Leo). Without names, the 48 ERA5-init projections don't pool with their ESM-init twins or with Goelzer 2025. |
| 4 | **High** | **ACCESS1.3's ERA5-init runs are labelled SSP5-8.5.** ACCESS1.3 is a CMIP5 model; the README codes "rcp85" and "ssp585" both as f85, and the ESM-init branch already reads RCP8.5. They should be RCP8.5. |
| 5 | Medium | **Control runs are in the app as a scenario.** "the SMB and ST anomalies are set to 0", holding the 1960–1989 reference climate from 2015. That's 33 runs, a quarter of the study. ISMIP6's control runs aren't loaded for any other study. Decision R2. |
| 6 | Medium | **The two initialization branches share one model row.** ERA5-init has one ERA5-forced spin-up and anomaly forcing. ESM-init has one spin-up per ESM and absolute ESM forcing. R25 is built around comparing the two. Decision R1. |
| 7 | Medium | Climate-model names aren't the official CMIP identifiers (e.g. "ACCESS1.3", "CESM2-CMIP6", "UKESM1-0-LL-CMIP6", "NorESM2"). Map them as Goelzer 2025 does (rules D5, G5). |
| 8 | Low | Ocean forcing, regional climate model and resolution aren't encoded. R25 follows "the ISMIP6 parameterization (Goelzer et al., 2020b)" for retreat, which is **Standard**. All ESM forcing is "dynamically downscaled over Greenland with the Modèle Atmosphérique Régional (MAR) v3.12", so **MARv3.12**. 4 km is **≤5 km**. |

**What was verified.**
- **Tier labels.** The κ tiers follow the physics. For the same branch, ESM and scenario, High loses more ice by 2100 than Medium in 31 of 32 pairs, and Medium more than Low in 32 of 32. The README codes `["med"]="50" ["high"]="25" ["low"]="75"`, as κ percentiles.
- **Run counts.** They match Table A1: 16 ESM–scenario projections per branch and tier.
- **Time stamps.** Each value is the state on 1 January of its year label, as in Coulon 2024.
- **Seam.** The historical run's last value (1 Jan 2015) and the projection's first (1 Jan 2016) join without a gap or duplicate.
- **Mass variable.** The loaded variable is absolute mass above flotation, not a change. That's harmless, because the app zeroes series at the window start and fits rates as slopes.
- **Resolution scope.** Only 4 km runs are loaded. ERA5-init also exists at 8 and 16 km, which are resolution-sensitivity copies of the same runs.

---

## 2. Experiments

- **Scenarios:** SSP1-2.6, SSP2-4.5, SSP5-8.5, and RCP8.5 for ACCESS1.3 (Table A1). "Control" holds the reference climate (R2).
- **Climate models.**
  - Ten ESMs: nine CMIP6 and one CMIP5.
  - Proposed CMIP identifiers, matching Goelzer 2025:

    | Archive name | Encoded |
    |---|---|
    | ACCESS1.3 | ACCESS1-3 |
    | CNRM-CM6 | CNRM-CM6-1 |
    | CNRM-ESM2 | CNRM-ESM2-1 |
    | CESM2-CMIP6 | CESM2 |
    | UKESM1-0-LL-CMIP6 | UKESM1-0-LL |
    | NorESM2 | NorESM2-MM |
    | CESM2-Leo | CESM2-Leo (kept distinct, G5) |
    | UKESM1-0-LL-Robin | UKESM1-0-LL-Robin (kept distinct, G5) |
    | MPI-ESM1-2-HR | unchanged |
    | IPSL-CM6A-LR | unchanged |

  - ERA5-init: the 1960–2014 historical run is forced by "the ERA5 reanalysis, downscaled with MAR". Projections use "anomalies with respect to the annual mean SMB and ST of the reference period (1960–1989) from the same ESM". The ESM is what forces the 2015–2023 window, so `climate_model` is that ESM, with "1960–2014: ERA5 (MAR)" in the hover.
  - ESM-init: forced "with absolute values of SMB and ST from the respective ESM" throughout.
- **Regional climate model:** MARv3.12 for every run.
- **Ocean forcing:** Standard. R25: "We follow Slater et al. (2019, 2020) and use the ISMIP6 parameterization".
- **Ocean sensitivity:** "three different values for κ covering the median, 25th percentile, and 75th percentile ... referred to as medium, high, and low sensitivity". Encoded **Medium / High / Low**.
- **Ice-shelf collapse:** none. "All floating ice is assumed to calve immediately".

---

## 3. Model configuration (CISM)

| Characteristic | Sources | Encoding |
|---|---|---|
| Ice model | "the Community Ice Sheet Model (CISM) (Lipscomb et al., 2019), which is a 3D thermomechanically coupled higher-order model" | CISM |
| Stress balance | "solving a depth-integrated viscosity approximation of the Stokes equations (Goldberg, 2011)" (DIVA) | **Higher-order** |
| Sliding law | See the notes below the table. | **Budd**, cascade step 2 (confidence medium; see the notes) |
| Initialization | 5000-year spin-up under 1960–1989 SMB; "basal friction parameters are calibrated to nudge the ice surface elevation towards observed conditions following Pollard and Deconto (2012)"; then "relax for an additional 1000 years, for the 4 km grid resolution"; "We assign the resulting ice sheet geometry to the beginning of 1960" | **Hybrid (thickness-nudged)** (D1), initial year 1960 |
| Resolution | "a horizontal grid resolution of 4, 8, and 16 km"; "All projections using ESM-init are run at 4 km grid resolution, while all ERA5-init experiments are run at three different grid sizes"; the app loads 4 km only | **≤5 km** (4 km) |

**Notes on the sliding law.**
- R25 Eq. 1: "a power law for basal sliding based on Weertman (1957): τb = k^(−1/p) N^(q/p) ub^(1/p) ... N is the effective pressure ... p = 3 and q = 1".
- Basal stress scales with N^(1/3). That makes it a Budd-type law (cascade step 2: friction depends on effective pressure), even though the paper calls it a Weertman power law.
- "Basal meltwater is removed immediately", so N is likely close to the ice overburden pressure. Friction still depends on N, and so on ice thickness.
- The friction coefficient is nudged; the law isn't Coulomb-limited.

---

## 4. Recommended encoding

**Model rows** (decision R1, Option A): `Rahlves2025 / CISM (ERA5-init)` and `Rahlves2025 / CISM (ESM-init)`. Both share these values:

| Field | Value |
|---|---|
| ice_model | CISM |
| stress_balance | Higher-order |
| resolution | ≤5 km (4 km) |
| initialization | Hybrid (thickness-nudged), 1960 |
| sliding_law | Budd (p = 3, q = 1) |

They differ only in their initialization descriptions: one ERA5 spin-up, or one spin-up per ESM.

**Experiments:**

| Field | Value |
|---|---|
| climate_model | CMIP identifier of the projection ESM |
| scenario | as in Table A1, ACCESS1.3 as RCP8.5 |
| regional_climate_model | MARv3.12 |
| ocean_forcing | Standard |
| ocean_sensitivity | Low / Medium / High from the tier |
| ice_shelf_fracture | No |
| hover | branch, historical forcing, κ percentile |

---

## 5. Decisions (applied 2026-11-23)

The user chose **Option A** for both decisions.
- **R1:** two model rows, `CISM (ERA5-init)` and `CISM (ESM-init)`, each with 48 runs (16 projections × 3 retreat sensitivities). Rahlves now counts as two configurations.
- **R2:** the 33 control runs are left out. The loader skips them, and `external_sources_rahlves2025_gis.csv.gz` was regenerated: 96 runs. The ESM-init values are unchanged.

**Also applied:**
- **Loader fixes:** the ERA5-init models are decoded from the archive README; ACCESS1.3 is RCP8.5 in both branches; model names map to CMIP identifiers as for Goelzer 2025.
- **Encodings:** Higher-order, Budd, Hybrid (thickness-nudged), ≤5 km, Standard, MARv3.12.
- **Files:** `Rahlves2025_models.csv` and `Rahlves2025_experiments.csv`.

What each decision weighed:

- **R1. Configurations.**
  - **Option A:** two model rows, ERA5-init and ESM-init. The branches differ in procedure (one ERA5 spin-up with anomaly forcing, against per-ESM spin-ups with absolute forcing), which is the paper's central comparison. Rahlves then counts as two configurations.
  - **Option B:** one row, as now.
  - **Option C:** eleven rows, one per spin-up (ERA5 plus each ESM).
  - Recommendation: **A**. Unlike Coulon's two initial states (C3), which were a sampled input within one procedure, these are two different procedures.
- **R2. Control runs.**
  - **Option A:** leave the 33 control runs out of the app. They hold the 1960–1989 climate and aren't climate projections, and no other study's control runs are loaded.
  - **Option B:** keep them as scenario "Control" but leave them out of the bias analysis.
  - **Option C:** keep them as now.
  - Recommendation: **A**.
