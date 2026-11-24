# Metadata review 5: DeConto et al. (2021), PSU3D-ICE Antarctic ensembles

**Paper:** DeConto, R. M., Pollard, D., Alley, R. B., et al. (2021). The Paris Climate Agreement and future sea-level rise from Antarctica. *Nature* 593, 83–89. https://doi.org/10.1038/s41586-021-03427-0. The main text, Methods, Extended Data captions and Supplementary Information were read. Citations to the paper below are "DC21".

**Data:** the paper's Source Data for Figure 1 (`41586_2021_3427_MOESM2_ESM.xlsx`).
- It has four sheets, "+1.5 ºC", "+2.0 ºC", "+3.0 ºC" and "RCP8.5".
- Each sheet has the same 109 runs, named `param_<a>_<b>`, as annual Antarctic sea-level contributions in metres, zeroed at 2000.
- The app uses 2000–2100: 436 runs.

**Rules applied:** those of reviews 1–4 (`encoding_rules.csv`).

**Supporting sources:**
- Pollard & DeConto (2012a), GMD 5, 1273, the model description: DC21 ref. 51.
- Pollard & DeConto (2012b), The Cryosphere 6, 953, the inverse method: DC21 ref. 52.
- DeConto & Pollard (2016), review 04: DC21 ref. 8.

---

## 1. Summary

| # | Severity | Finding |
|---|---|---|
| 1 | **High** | The model metadata (Weertman sliding, Spin-up) was entered without a source. DC21 gives Power law (Weertman) sliding, Hybrid (thickness-nudged) initialization, Shallow stress balance and 10 km resolution (§3). |
| 2 | Medium | `climate_model` reads "Not reported in source data". DC21 reports it. The ocean comes from transient NCAR CCSM4 runs (RCP4.5 or RCP8.5). The atmosphere comes from **RegCM3** snapshots at fixed CO2 levels, whose sea-surface temperatures and sea ice come from the same CCSM4 runs (§2). Decision DC1. |
| 3 | Medium | **All 109 members were selected partly by matching IMBIE.** The first test is "the average annual mass change dM/dt from 1992–2017 ... provided by the IMBIE assessment"; 33 of 196 fail it. The app compares rates with IMBIE3 over 2015–2023, which overlaps the training period for three years. Decision DC3. |
| 4 | Low | The hover labels the parameters "CREVLIQ, CLIFVMAX", a name borrowed from another dataset. DC21 calls them CALVLIQ (also CREVLIQ) and VCLIFF. Rename to the paper's names. |
| 5 | Low | The +1.5 °C and +2.0 °C scenarios aren't RCPs. Both follow RCP4.5 until warming reaches the target (2040, 2060) and then hold the atmosphere and ocean fixed. +3.0 °C follows RCP8.5 until 2070, then holds. This should be shown in the hover (§2). |

**What was verified.**
- **Parameter reading.** `param_<a>_<b>` is a = CALVLIQ, b = VCLIFF. The 109 runs take CALVLIQ 15–195 (step 15) and VCLIFF 2–13 km/yr.
  - Their means are 107.1 and 7.65 and their medians 105 and 7. Extended Data Table 1 gives "average calibrated parameter values ... CALVLIQ = 107 m−1 yr2 and VCLIFF = 7.7 km yr−1. The corresponding median values are 105 m−1 yr2 and 7 km yr−1".
- **Ensemble.** The 109 runs are the calibrated ensemble: 196 runs, then 163 after IMBIE, 119 after the LIG and "n = 109" after the Pliocene test.
- **2100 values.** The 2100 median and 17th–83rd percentiles of every sheet reproduce Table 1 exactly:
  - +1.5 °C: 0.08 (0.06–0.10) m
  - +2.0 °C: 0.09 (0.07–0.11) m
  - +3.0 °C: 0.15 (0.08–0.27) m
  - RCP8.5: 0.34 (0.20–0.53) m
- **Sheet differences.** The +1.5 °C and +2.0 °C sheets are identical until 2036, as expected since both follow RCP4.5.
  - The +3.0 °C and RCP8.5 sheets differ by a few mm from 2001 (≤ 3.8 mm by 2023), although DC21 describes +3.0 °C as RCP8.5-based until 2070.
  - This doesn't affect any encoding, and the data are used as published.

---

## 2. Experiments

- **Scenarios** (DC21 Methods, "Future ensembles"):
  - **+1.5 °C:** "initially follow a RCP4.5 emission trajectory ... evolves freely until 2040 ... the atmosphere and ocean forcings are fixed (maintained) at their 2040 levels".
  - **+2.0 °C:** "also based on RCP4.5, but warming is allowed to evolve until 2060".
  - **+3.0 °C:** "our +3 °C scenario (roughly representing the NDCs) is based on RCP8.5, with atmospheric and oceanic forcing fixed beyond 2070". Table 1 calls it "+3.0 °C (NDCs) ensemble".
  - **RCP8.5:** extended RCP8.5.
  - Encoding: keep the scenario labels (+1.5 °C, +2.0 °C, +3.0 °C, RCP8.5). The base pathway and fixing year go in the hover. RCP8.5 joins the "composite 8.5" grouping; the warming-level scenarios stay separate.
- **Period:** "Future ice-sheet simulations begin in 1950 to allow comparisons with observations over the satellite era" (DC21 Methods).
- **Ensemble.** "Our primary perturbed physics ensembles use a 14 × 14 matrix (n = 196) of CREVLIQ and VCLIFF". The ensemble is scored on three tests:
  - **IMBIE:** dM/dt 1992–2017, "equivalent to a GMSL change of 0.15–0.46 mm yr−1".
  - **LIG:** 3.1–6.1 m.
  - **Pliocene:** 11–21 m.

  The runs are unweighted ("weighting of individual model outcomes is avoided").
- **Atmospheric forcing.**
  - "RCM snapshots are run at 1950 and with increasing levels of effective CO2 (2, 4 and 8 times the preindustrial level) ... time-interpolated and log-weighted to match transient CO2 concentrations".
  - The RCM is RegCM3: ref. 79, Pal et al. 2007, "the ITCP RegCM3".
  - Unlike DP16, "sea surface temperatures and sea ice boundary conditions in the nested RCM come from the same transient NCAR CCSM4 runs that provide the time-evolving sub-surface ocean temperatures".
  - A uniform 2.9 °C warm correction is applied for the RCM's summer cold bias.
  - DC21 doesn't name the GCM supplying the RCM's lateral atmospheric boundaries; DP16 used GENESIS v3.
- **Ocean forcing.**
  - CCSM4, ref. 80: Shields & Kiehl 2016, "the half degree Community Climate System Model, CCSM4".
  - Temperatures at 400 m drive a quadratic melt relation (Pollard & DeConto 2012a), with "a single, uniform coefficient in the ocean melt relation based on a 625-member ensemble of simulations of WAIS retreat through the last deglaciation".
  - "A 1.5 °C sub-surface ocean temperature adjustment is used in the Amundsen Sea sector".
  - Encoding: the group's own parameterization, so **Open**, as in DP16 (decision DP3).
- **Ocean sensitivity.** There is one melt coefficient for all runs, calibrated by the group. Decision DC2.
- **Ice-shelf collapse:** none imposed. Hydrofracturing (CALVLIQ) and cliff calving (VCLIFF) are model physics, so `ice_shelf_fracture` = No.

---

## 3. Model configuration (PSU3D-ICE)

| Characteristic | Sources | Encoding |
|---|---|---|
| Ice model | "hybrid ice dynamics" (ref. 51, Pollard & DeConto 2012a), "Several advances relative to previous versions of the model" (DC21 Methods) | PSU3D-ICE |
| Stress balance | "The ice sheet–shelf model uses hybrid ice dynamics"; SSA equations with SIA (Pollard & DeConto 2012a: "a hybrid combination of the scaled shallow ice and shallow shelf approximations") | **Shallow** |
| Sliding law | See the notes below the table. | **Power law (Weertman)**, m = 2 |
| Initialization | See the notes below the table. | **Hybrid (thickness-nudged)** (D1 rule), initial year 1950 |
| Resolution | "A resolution of 10 km is used for continental simulations used in our main ensembles (Figs. 1–3). A nested 1-km grid is used for a select simulation over West Antarctica" | 10 km → **>5–10 km**. The 1-km nest isn't in the Figure 1 ensembles. |

**Notes on the sliding law.**
- DC21: "The model uses a standard Weertman-type basal sliding law [ref. 51], with basal sliding coefficients determined by an inverse method".
- Ref. 51 (Pollard & DeConto 2012a, Eqs. 10–11): u_b = C′ |τ_b|^(m−1) τ_b with "m = 2", ramped by basal temperature, and no effective pressure.
- Cascade step 4 applies.

**Notes on initialization.**
- "initial ice-sheet conditions (ice thickness, bed elevation, velocity, basal sliding coefficients and internal ice and bed temperatures) follow the same procedure as in ref. 8 and are identical in all simulations. Initialization involves a 100,000-kyr [sic] spin-up using observed mean annual ocean climatology and standard SeaRISE atmospheric temperature and precipitation fields."
- The sliding coefficients come from an inverse method "that iteratively matches model ice-surface elevations to observations".
- "Future ice-sheet simulations begin in 1950".

---

## 4. Recommended encoding

**Model row `DeConto2021 / PSU3D-ICE`:**

| Field | Value |
|---|---|
| review_status | reviewed |
| stress_balance | Shallow |
| resolution | >5–10 km (res_min_km = res_max_km = 10) |
| initialization | Hybrid (thickness-nudged), initial year 1950 |
| sliding_law | Power law (Weertman), m = 2, no effective pressure, high confidence |

**Experiments:**

| Field | Value |
|---|---|
| scenario | unchanged |
| climate_model | CCSM4 (ocean); RegCM3 CO2 snapshots (atmosphere) (DC1) |
| regional_climate_model | RegCM3 |
| ocean_forcing | Open |
| ocean_sensitivity | Medium (DC2) |
| ice_shelf_fracture | No |
| hover | CALVLIQ, VCLIFF (km/yr), base pathway, calibration (DC3) |

---

## 5. Decisions (applied 2026-11-23)

The user chose **Option A** for each decision.
- **DC1:** climate_model is **CCSM4 (ocean); RegCM3 CO2 snapshots (atmosphere)**. The loader (`utilities/external_sources.py`) writes the same label, and the bundled `external_sources_deconto2021_ais.csv.gz` was regenerated with it; its mass values are unchanged.
- **DC2:** ocean_sensitivity is **Medium**. `encoding_rules.csv` now names DC21 under Medium.
- **DC3:** DC21 stays in the bias analysis. Every run's hover states that members were selected partly by matching IMBIE 1992–2017.

The hover also gives each run's CALVLIQ, VCLIFF (renamed from "CREVLIQ, CLIFVMAX") and emissions pathway. The metadata is in `dash_app/metadata/DeConto2021_models.csv` and `DeConto2021_experiments.csv`.

What each decision weighed:

- **DC1. Climate model label.**
  - **Option A:** "CCSM4 (ocean); RegCM3 CO2 snapshots (atmosphere)". This keeps DC21 apart from ISMIP6 runs forced entirely by CCSM4 output, as DP4 did for DP16. It doesn't name a parent atmospheric GCM, because DC21 doesn't either.
  - **Option B:** "CCSM4". CCSM4 provides the ocean and the RCM's sea surface, and RegCM3 is already recorded as the regional climate model.
  - Recommendation: **A**.
- **DC2. Ocean sensitivity.**
  - **Option A:** **Medium**. The existing rule gives Medium to "the group's own calibration (open)", as for ISMIP6 open-framework runs. DC21 uses one group-calibrated coefficient for every run.
  - **Option B:** "Study-specific (see experiment)", as for DP16. DP16, however, varies its coefficient across runs (OCFAC 0.1–10×), and DC21 doesn't.
  - Recommendation: **A**.
- **DC3. Calibration against IMBIE.**
  - **Option A:** keep DC21 in the bias analysis. Record "members selected to match IMBIE 1992–2017 dM/dt (0.15–0.46 mm/yr GMSL)" in the hover and in this report.
  - **Option B:** leave DC21 out of the bias analysis, like the ice-shelf collapse runs.

  Recommendation: **A**, because:
  - The overlap with the app's 2015–2023 window is three years.
  - The target is a broad 25-year mean.
  - Excluding studies tuned to observations would also call ISMIP6's velocity-inversion models into question.
