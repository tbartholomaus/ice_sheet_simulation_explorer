# Metadata review 8: Aschwanden et al. (2019), PISM Greenland large ensemble

**Paper:** Aschwanden, A., Fahnestock, M. A., Truffer, M., Brinkerhoff, D. J., Hock, R., Khroulev, C., Mottram, R. and Khan, S. A. (2019). Contribution of the Greenland Ice Sheet to sea level over the next millennium. *Science Advances* 5, eaav9396. https://doi.org/10.1126/sciadv.aav9396 (open access; full text read from PMC6584365). Citations to the paper below are "A19".

**Archive:** Aschwanden, A. and Brinkerhoff, D. J. (2022), "Calibrated mass loss projections from the Greenland Ice Sheet", Arctic Data Center, https://doi.org/10.18739/A2KW57K4R (the data set of Aschwanden & Brinkerhoff 2022, GRL 49, e2022GL099058; "A&B22").
- **Runs:** `pism_scalars/aschwanden_et_al_2019_les_2008_norm.csv.gz`. The archive's metadata record was read; the GRL paper itself refused automated download (HTTP 403).
- **Parameter table:** `lhs_samples_gcm.csv` (A19's archive, doi:10.18739/A2222R58F, and the authors' code repository `pism/pism-gris`, `latin_hypercube/`).
- **Code that ran the ensemble:** `pism/pism-gris/stability/lhs_ensemble.py` and `latin_hypercube/draw_samples.py`.
- **A19's ensemble methods document**, Table 1 (cn.dataone.org urn:uuid:e9780902-…).

**What the app uses:** 1458 runs: 489 (RCP2.6), 491 (RCP4.5) and 478 (RCP8.5) of the 500 sampled members. Annual values from 2008, zeroed at 2008, trimmed at 2100.

**Rules applied:** those of reviews 1–6 (`encoding_rules.csv`).

---

## 1. Summary

| # | Severity | Finding |
|---|---|---|
| 1 | **High** | **The loader's provenance note is wrong about what the runs are.** It says A&B22 "draws large ensembles from" a neural-network emulator and that the loaded "les" file is such a draw. The runs are in fact real PISM simulations from A19 (note below the table). The study key, label ("Aschwanden 2019") and README citation are therefore right; the loader comment and docstring need correcting. Decision AS3. |
| 2 | **High** | **The model metadata was entered without a source**, and the initialization is wrong. The encoded "Spin-up" is wrong: A19 starts from a paleo spin-up and then a 100-year relaxation that keeps "the ice surface close to observations using a flux correction". That is **Hybrid (thickness-nudged)** (rule D1), exactly as the same initial state is encoded for UAF-PISM in Goelzer 2020. "Pseudo-plastic" isn't a value in the vocabulary: A19's law is pseudo-plastic with a Mohr–Coulomb yield stress, so **Coulomb-limited**. Stress balance and resolution weren't encoded: hybrid SIA+SSA (**Shallow**) at **1.8 km** (**≤5 km**). |
| 3 | Medium | **The regional climate model isn't encoded.** A19 forces PISM with a fixed 2000–2015 climatology from HIRHAM5 (driven by ERA-Interim, not by a GCM), shifted each year by one **Greenland-mean** CMIP5 temperature anomaly. Decision AS1. |
| 4 | Medium | **Ocean forcing and sensitivity aren't encoded.** A19 uses its own submarine-melt parameterization driven by air-temperature anomalies. Each member draws a low, mid or high regime for the present-day spatial pattern (OCM) and, independently, for the temporal response to warming (OCS). Decision AS2. |
| 5 | Low | The protocol hover says only "Uncalibrated (prior ensemble)". It should say this is A19's ensemble, unweighted: A&B22's calibration, which uses mass change 2010–2020, isn't applied. |
| 6 | Low | 22, 9 and 11 of the 500 members are missing for RCP8.5, 4.5 and 2.6 in the archive (presumably runs that failed). The data are used as published. |

**What finding 1 rests on.**
- The archive's metadata record says the scalar files are "prepared with the Parallel Ice Sheet Model (PISM) ... for both the ensemble using the Prior and the Posterior (calibrated) parameter distribution". The emulator (trained on ~1000 surface-speed runs) is used only to calibrate the parameters.
- The "les" file is A19's large ensemble:
  - its file name;
  - its 2008–3007 span;
  - its single resolution, 1800 m;
  - its member ids 0–499, which match A19's `lhs_samples_gcm.csv`;
  - its 16th–84th percentiles, which reproduce A19 Table 1 (§1, "What was verified").

**What was verified.**
- **The runs are A19's ensemble.** 16th–84th percentile GMSL relative to 2008, in cm (computed at 362.5 Gt/mm; A19 Table 1 in brackets):

| Year | RCP2.6 | RCP4.5 | RCP8.5 |
|---|---|---|---|
| 2100 | 6–19 (5–19) | 8–23 (8–23) | 14–34 (14–33) |
| 2300 | 17–55 (17–55) | 36–98 (35–97) | 94–377 (94–374) |
| 3000 | 59–192 (59–191) | 198–419 (197–416) | 536–734 (538–728) |

- **The GCM code order.** The loader's code-to-name map ({0: GISS-E2-H, 1: GISS-E2-R, 2: IPSL-CM5A-LR, 3: MPI-ESM-LR}) was an inference from A19's prose. It is now confirmed by the authors' own run script (`stability/lhs_ensemble.py`), which uses `gcm_dict = {0.0: "GISS-E2-H", 1.0: "GISS-E2-R", 2.0: "IPSL-CM5A-LR", 3.0: "MPI-ESM-LR"}` to choose each member's forcing file. As a sanity check, code 2 (IPSL-CM5A-LR, the warmest of the four) gives the largest RCP8.5 losses (median 29 cm by 2100 against 16–22 cm for the others).
- **The rest of the parameter table** matches A19's Table 1 design:
  - GCM is uniform over 4 models (125 members each).
  - OCM, OCS and TCT each take −1/0/1 (low/mid/high) about 167 times.
  - PPQ (q) spans 0.25–0.95.
- **Bundled data:** 1458 runs, 2008–2100, zeroed at 2008, mass loss negative. No joins or seams: each run is one continuous PISM run from 2008.

---

## 2. Experiments

- **Scenarios:** RCP2.6, RCP4.5 and RCP8.5, for each of 500 Latin-hypercube members ("for each RCP scenario, we produced 500 parameter combinations"). Correct.
- **Climate forcing.**
  - A19: "we adjust T0(x, y, t) by scalar temperature anomalies ΔTair(t) derived from a GCM ... we select the four GCMs that extend to year 2300: GISS-E2-H, GISS-E2-R, IPSL-CM5A-LR, and MPI-ESM-LR ... we first extract monthly near-surface temperature averaged over the domain containing Greenland ... subtract the 2006–2015 monthly means".
  - T0 and P0 are HIRHAM5 2000–2015 means ("forced at the lateral boundaries using ... ERA-Interim"); precipitation scales by ω = 5–7 % per K.
  - Encoding: `climate_model` is the member's CMIP5 GCM (official identifiers, already correct).
- **Regional climate model:** HIRHAM5 supplies only the fixed present-day climatology; no GCM is downscaled. Decision AS1.
- **Ocean forcing.**
  - A19's own parameterization of submarine (frontal and sub-shelf) melt. A present-day rate varies with latitude (OCM low/mid/high: (ṁ71, ṁ80) = (300, 10), (400, 20) or (500, 30) m/yr). It scales with warming as 1 + C(aΔTair)^α ΔTair^β (OCS low/mid/high: (α, β) = (0.5, 1.0), (0.54, 1.17) or (0.85, 1.61)).
  - A19 "assume that ocean thermal forcing Th increases by 1 K for every 1-K increase in air temperature".
  - Encoding: the group's own scheme, so **Open**.
- **Ocean sensitivity:** two independent tiers per member (OCM, OCS) on A19's own scale. Decision AS2.
- **Calving:** a von Mises stress criterion (σmax, VCM) plus a minimum floating-ice thickness (h_min, TCT low/mid/high). This is model physics, not an imposed collapse, so `ice_shelf_fracture` = No.
- **Calibration:** none in A19. A&B22 later calibrated against surface speeds and against "cumulative mass change from 2010 until 2020", but the app uses the uncalibrated ("prior") file, so no member is selected or weighted on observations.

---

## 3. Model configuration (PISM)

| Characteristic | Sources | Encoding |
|---|---|---|
| Ice model | "the open-source PISM ... thermomechanically coupled and polythermal" (A19 Methods) | PISM |
| Stress balance | "uses a hybrid stress balance" (A19); `--stress_balance ssa+sia` default in `lhs_ensemble.py` | **Shallow** |
| Sliding law | See the notes below the table. | **Coulomb-limited**; q sampled 0.25–0.95 (median 0.6); effective pressure yes (till) |
| Initialization | See the notes below the table. | **Hybrid (thickness-nudged)** (D1), initial year 2008 |
| Resolution | "the ensemble simulations using resolutions of ... 1800 m" (A19); archive column `Resolution (m)` = 1800 for every row | res_min = res_max = 1.8 km → **≤5 km** |

**Notes on the sliding law.**
- A19: "a pseudoplastic power law relates bed-parallel shear stress and sliding", τb = −τc ub/(\|ub\|^(1−q) u0^q).
- "We assume that yield stress τc is proportional to effective pressure N ['Mohr-Coulomb criterion'] τc = (tan ϕ) N".
- The exponent q is an ensemble parameter: PPQ, a truncated normal with μ = 0.6 and σ = 0.2.
- Cascade step 1 applies: the basal stress is capped by a yield stress. This is the same encoding as UAF-PISM in Goelzer 2020.

**Notes on initialization.**
- A19: "Ice sheet initial conditions are provided by the 'calibrated' experiments in (12) [Aschwanden et al. 2016], which includes a 125-ka paleoclimate simulation followed by calibration of ice dynamical parameters to minimized differences between observed and simulated surface velocities".
- Then: "we then ran a 100-year long relaxation simulation ... but we kept the ice surface close to observations using a flux correction. The result was an initial state that is both close to the observed geometry and surface speeds of 2008."
- The calibration tunes a few scalar parameters (q, E, the till friction angle against bed elevation), not a spatial inversion of basal friction, so this isn't "Data assimilation (velocity inversion)".
- Goelzer 2020's UAF-PISM rows use the same initial state and are encoded Hybrid (CYC/NDs).

---

## 4. Recommended encoding

**Model row `Aschwanden2019 / PISM`:**

| Field | Value |
|---|---|
| review_status | reviewed |
| ice_model | PISM |
| stress_balance | Shallow (code Hybrid) |
| resolution | ≤5 km (res_min_km = res_max_km = 1.8) |
| initialization | Hybrid (thickness-nudged), initial year 2008 |
| sliding_law | Coulomb-limited, "q = 0.25–0.95 (sampled, median 0.6)", effective pressure yes (till), high confidence |

**Experiments file:** transcribe all 11 Latin-hypercube parameters per member from `lhs_samples_gcm.csv`, decoded per `lhs_ensemble.py`:

| Parameter | Meaning |
|---|---|
| GCM | 0–3 → model name |
| PRS | ω |
| FICE | ice melt factor |
| FSNOW | snow melt factor |
| RFR | refreezing proportion |
| OCM | −1/0/1 → low/mid/high |
| OCS | −1/0/1 → low/mid/high |
| TCT | −1/0/1 → low/mid/high |
| VCM | σmax, MPa |
| PPQ | q |
| SIAE | E |

**Experiment fields:**

| Field | Value |
|---|---|
| climate_model | the member's CMIP5 GCM (unchanged) |
| scenario | RCP2.6, RCP4.5, RCP8.5 (unchanged) |
| regional_climate_model | None (AOGCM output used directly) (AS1) |
| ocean_forcing | **Open** |
| ocean_sensitivity | Study-specific (see experiment) (AS2) |
| ice_shelf_fracture | No |
| hover | OCM and OCS tiers, q, the forcing description (AS1) and the protocol "Aschwanden et al. 2019 ensemble, unweighted (A&B22's calibration not applied)". The other parameters are "(file)" only, to keep the hover short. |

---

## 5. Decisions (applied 2026-11-23)

The user chose the recommended option for each decision, and approved the encodings in §3–4.
- **AS1:** regional_climate_model is **None (AOGCM output used directly)**. The hover's "Climate forcing" line reads "HIRHAM5 climatology + Greenland-mean GCM ΔT".
- **AS2:** ocean_sensitivity is **Study-specific (see experiment)**. The hover shows each run's OCM and OCS tiers and q.
- **AS3:** the comment and docstring of `load_aschwanden2022_gis()` in `utilities/external_sources.py` are corrected: real PISM runs of A19's ensemble, and the GCM code order is confirmed by `pism-gris/stability/lhs_ensemble.py`. No code or identifiers changed.

The model row is reviewed: Shallow, Coulomb-limited (q 0.25–0.95), Hybrid (thickness-nudged) from 2008, ≤5 km (1.8 km), ocean forcing Open, no ice-shelf collapse. `Aschwanden2019_experiments.csv` transcribes all 11 Latin-hypercube parameters for each of the 1458 runs. OCM, OCS and q are shown in the hover, and the others are kept in the file only. The metadata is in `dash_app/metadata/Aschwanden2019_models.csv` and `Aschwanden2019_experiments.csv`.

What each decision weighed:

- **AS1. Regional climate model.**
  - **Option A:** "None (AOGCM output used directly)", as for ISMIP6-AIS and Coulon 2024 (GCM anomalies added to a present-day climatology). The hover would add "HIRHAM5 2000–2015 climatology shifted by a Greenland-mean GCM temperature anomaly".
  - **Option B:** "HIRHAM5". That would group A19 with runs whose GCM forcing HIRHAM5 actually downscaled, which isn't the case here.
  - Recommendation: **A**.
- **AS2. Ocean sensitivity.**
  - **Option A:** "Study-specific (see experiment)", with the OCM and OCS tiers in the hover. A19's low/mid/high regimes are its own scale on two independent dimensions, not percentiles of the ISMIP6 κ distribution. This is the treatment DeConto & Pollard 2016 and Coulon 2024 got.
  - **Option B:** map OCS (the response to warming, the closest analogue of κ) to Low/Medium/High, so A19 pools with the ISMIP6-GIS and PROTECT tiers. OCM would be shown in the hover.
  - Recommendation: **A**.
- **AS3. Loader provenance.**
  - **Option A:** correct the comment and docstring of `load_aschwanden2022_gis()`. The "les" file holds real PISM runs of A19's ensemble, as archived with A&B22, and the GCM code order is confirmed by `pism-gris/stability/lhs_ensemble.py`. Function and file names stay as they are.
  - **Option B:** additionally rename the "2022" identifiers (`load_aschwanden2022_gis`, `external_sources_aschwanden2022_gis.csv.gz`, `aschwanden2022_gis` in `data.py`) to "2019".
  - Recommendation: **A**. The names don't reach the app's users, and renaming touches the data file and the notebook.
