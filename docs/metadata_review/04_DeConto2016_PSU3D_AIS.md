# Metadata review 4: DeConto & Pollard (2016), PSU3D-ICE Antarctic ensemble

**Paper:** DeConto, R. M. and Pollard, D. (2016). Contribution of Antarctica to past and future sea-level rise. *Nature* 531, 591–597. https://doi.org/10.1038/nature17145. The full text was read, including the Methods and the Extended Data captions. Citations to the paper below are "DP16".

**Data:** Kopp et al. (2017), *Earth's Future* 5, 1217–1233, https://doi.org/10.1002/2017EF000663, Data Set S1. It is bundled unmodified as `eft2271-sup-0002-2017ef000663-ds01.tsv` and converted when the app starts by `dash_app/data.py`.
- Before this review, the app used the AIS-total rows "RCP 2.6", "RCP 4.5" and "RCP 8.5": 29 members × 3 scenarios = 87 runs. Decadal values are interpolated to annual.
- Data Set S1 also has "RCP 2.6 (BC)", "RCP 4.5 (BC)" and "RCP 8.5 (BC)" rows. They have the same 29 members, with the same OCFAC, CREVLIQ, VCLIF, Pliocene and LIG values, 29 per scenario.
- Since decision DP1, the app uses the "(BC)" rows.

**Rules applied:** those of reviews 1–3 (`encoding_rules.csv`).

**Supporting sources:**
- Pollard & DeConto (2012a), GMD 5, 1273, the model description: DP16 ref. 27.
- Pollard & DeConto (2012b), The Cryosphere 6, 953, the inverse method: DP16 ref. 51.
- DeConto et al. (2021), Nature 593, 83, which uses the same model.

---

## 1. Summary

| # | Severity | Finding |
|---|---|---|
| 1 | **High** | **The app uses the Data Set S1 rows *without* DP16's ocean-temperature correction, so it isn't showing DP16 as published.** See the note below the table. Decision DP1. |
| 2 | **High** | The model metadata (Weertman sliding, Spin-up) was entered without a source. DP16 and the papers it cites give Power law (Weertman) sliding, Hybrid (thickness-nudged) initialization, Shallow stress balance and 10 km resolution (§3). |
| 3 | Medium | `climate_model` reads "Not reported in Data Set S1". DP16 has **two different climate models**: the ocean comes from CCSM4 RCP runs; the atmosphere from GENESIS v3 with the RegCM3 regional model. The atmospheric forcing is **not** RCP output from a CMIP GCM (§2). Decision DP4. |
| 4 | Medium | Ocean forcing and ocean sensitivity have no ISMIP6 counterpart. DP16 uses its own quadratic melt parameterization; OCFAC is a multiplier on its coefficient. Decision DP3. |
| 5 | Low | The 29 members are DP16's *second* calibration, with Pliocene targets of 5–15 m. DP16's headline numbers come from the 15-member set with targets of 10–20 m. Both are DP16 ensembles, and Kopp et al. 2017 chose the 29. Record it in the protocol hover text; no change to the runs (§2). |
| 6 | Low | The 29 members differ only in three physical parameters (OCFAC, CREVLIQ, VCLIF), which are already recorded per experiment. Stress balance, sliding, initialization and resolution are the same for all, so one model row describes them all. |

**Note on finding 1.**
- **What DP16 did:** "We find the addition of 3 °C to the CCSM4 ocean temperatures increases melt rates to 25–30 m yr−1 ... this substantially improves grounding-line positions in the Amundsen Sea (Pine Island Glacier in particular) from 1950 to 2015" (DP16 Methods).
- **DP16's headline results include the correction.** "Adding the ocean temperature correction ... further increases the 2100 projections in RCP2.6, RCP4.5 and RCP8.5 to 16 ± 16 cm, 58 ± 28 cm and 114 ± 36 cm". Both ensembles were run, "with and without the bias correction".
- **Kopp et al. 2017 use the corrected runs as DP16.** They show the uncorrected runs only as a sensitivity case.
- **Which rows are which.** Kopp's own numbers identify the rows unambiguously:
  - Plain rows, 2050: RCP8.5 median −28 mm (range −91 to +122) and RCP2.6 median −19 mm (−99 to +62). These reproduce Kopp's "Without this correction ... RCP 8.5: median −3 cm, range of −9 to +12 cm, RCP 2.6: median −2 cm, range of −10 to +6 cm".
  - "(BC)" rows, 2050: RCP8.5 median +52 mm. This reproduces "a median contribution to GMSL of about +5 cm under DP16".
- **Why the plain rows were chosen.** The 2026-10-08 choice rests on a loader note that calls the correction "Kopp et al.'s". It is DP16's own.

---

## 2. Experiments

- **Scenarios:** RCP2.6, RCP4.5 and RCP8.5, extended to 2500 (DP16). Correct.
- **Period:**
  - "The simulations begin in 1950 to provide some hindcast spinup, and are run for 550 years to 2500" (DP16).
  - Data Set S1 starts at 2000 = 0, and the app shifts it to its 2015 baseline.
- **Ensemble.**
  - Three parameters, four values each, give 64 combinations:
    - **OCFAC:** the sub-ice-shelf melt coefficient, "0.1, 1, 3 and 10 times the default value of 0.224 m yr−2 °C−2".
    - **CREVLIQ:** hydrofracturing, "0 m, 50 m, 100 m and 150 m per (m yr−1)−2".
    - **VCLIF:** maximum cliff-failure rate, "0 km yr−1, 1 km yr−1, 3 km yr−1 and 5 km yr−1".
  - The 64 were run for the Pliocene and the LIG.
  - With targets of 10–20 m (Pliocene) and 3.6–7.4 m (LIG), "15 out of 64" pass; these give DP16's headline numbers.
  - With "a lower accepted range for Pliocene sea-level rise (5–15 m) ... (29 out of 64 passed this test)".
  - Data Set S1, and therefore the app, uses the **29-member** set. The current protocol text ("passing Pliocene (5-15 m) and LIG tests") is correct. It should add that this is DP16's alternative, lower-target calibration.
- **Atmospheric forcing.**
  - It is not taken from an RCP GCM run. "The nested GCM–RCM is run to equilibrium with 1 × PAL, 2 × PAL, 4 × PAL and 8 × PAL CO2 ... the climate at any time is the average of the two appropriate surrounding RCM solutions, weighted according to the logarithm of the concentration of CO2" (DP16 Methods).
  - The GCM is GENESIS v3 and the RCM is RegCM3: "decadal averages of meteorological fields from the RegCM3 RCM ... The RCM uses a 40-km grid ... nested within the GENESIS v3 Global Climate Model".
  - The forcing is anomaly-corrected against ALBMAP.
- **Ocean forcing.**
  - "Ocean temperatures in the RCP scenarios are provided by high-resolution (0.5° atmosphere and 1° ocean) NCAR CCSM4 ocean model output ... run to 2300", held at 2300 values afterwards.
  - Temperatures at 400 m drive DP16's quadratic melt parameterization (Pollard & DeConto 2012a, eq. 17), reduced in confined embayments.
  - This is not an ISMIP6 protocol. Decision DP3.
- **Ocean sensitivity:** OCFAC is a multiplier on DP16's melt coefficient. It is not a percentile of the ISMIP6 γ0 or κ distributions. Decision DP3.
- **Ice-shelf collapse:** none imposed. Hydrofracturing and cliff failure are model physics, so `ice_shelf_fracture` = No (they're recorded per run as CREVLIQ and VCLIF).

---

## 3. Model configuration (PSU3D-ICE)

| Characteristic | Sources | Encoding |
|---|---|---|
| Ice model | "an established ice-sheet model, with hybrid ice dynamics following the formulation described in ref. 27" (DP16; ref. 27 = Pollard & DeConto 2012a, PSU3D-ICE) | PSU3D-ICE |
| Stress balance | "hybrid ice dynamics" (DP16). Pollard & DeConto 2012a: "a hybrid combination of the scaled shallow ice and shallow shelf approximations" | **Shallow** |
| Sliding law | See the notes below the table. | **Power law (Weertman)**, m = 2 (u ∝ τ²) |
| Initialization | See the notes below the table. | **Hybrid (thickness-nudged)**: a spin-up nudged to observed thickness and surface elevation (D1 rule). Initial year 1950. |
| Resolution | "A grid resolution of 10 km is used for all simulations" (DP16 Methods) | res_min = res_max = 10 km → **>5–10 km** (10 km is a bin edge and belongs to the lower bin) |

**Notes on the sliding law.**
- **DP16 doesn't restate the law.** It says only "Basal sliding coefficients are determined by an inverse method".
- **Ref. 27 (Pollard & DeConto 2012a), Sect. 2.4.** Basal sliding is "treated ... by a standard drag law", u_b = C′ |τ_b|^(m−1) τ_b with m = 2.
  - C′ ramps from no sliding below a homologous basal temperature of −3 °C to full sliding at the melt point.
  - C(x,y) comes from the inverse method.
  - There is no effective-pressure term.
- **Other source.** DeConto et al. 2021: "standard Weertman-type basal sliding law".
- **Cascade.** Not Coulomb-limited, no N dependence, m = 2 ≠ 1, so step 4 applies.

**Notes on initialization.**
- **DP16 Methods quotes:**
  - "Ice-sheet initial conditions and basal sliding coefficients are provided by a 100-kyr inverse simulation following the methodology in ref. 51, using mass-balance forcing provided by a bias-corrected RCM climatology and modern observed ocean temperatures. In the inverse procedure, basal sliding coefficients under modern grounded ice are adjusted iteratively to reduce the misfit with observed ice thickness, with grounding-line positions fixed to observed locations."
  - Also: "iteratively matching ice-surface elevations to observations until a quasi-equilibrium is reached".
- **Start of the run:** "The simulations begin in 1950".

---

## 4. Recommended encoding

**Model row `DeConto2016 / PSU3D-ICE`:**

| Field | Value |
|---|---|
| review_status | reviewed |
| stress_balance | Shallow |
| resolution | >5–10 km (res_min_km = res_max_km = 10) |
| initialization | Hybrid (thickness-nudged), initial year 1950 |
| sliding_law | Power law (Weertman); m = 2; no effective pressure; confidence high (stated in the model description the paper cites) |

**Experiments:**

| Field | Value |
|---|---|
| scenario | RCP2.6, RCP4.5, RCP8.5 (unchanged) |
| climate_model | CCSM4 (ocean); GENESIS v3 (atmosphere) (DP4) |
| regional_climate_model | RegCM3 |
| ocean_forcing | Open (DP3) |
| ocean_sensitivity | Study-specific (see experiment) (DP3) |
| ice_shelf_fracture | No |
| hover | OCFAC, CREVLIQ, VCLIF, Pliocene and LIG values, and the protocol; adds "ocean correction: +3 °C Amundsen–Bellingshausen" if DP1 switches rows |

---

## 5. Decisions (applied 2026-11-23)

The user chose the recommended option of each decision.
- **DP1 → A:** the app and notebook loaders (`dash_app/data.py`, `utilities/external_sources.py`) now read the "(BC)" rows. Experiment names are unchanged (`RCP85_m01` ...), and the experiments file's verbatim Scenario column ("RCP 8.5 (BC)") records the correction.
- **DP3:**
  - Ocean forcing is **Open**.
  - Ocean sensitivity is a new value, **Study-specific (see experiment)**, added to `encoding_rules.csv` after PIGL.
  - The hover shows each run's OCFAC, CREVLIQ, VCLIF, Pliocene and LIG values.
- **Later change (user request, 2026-11-23):** the Protocol line, which named the 29-member calibration and the +3 °C ocean correction, is kept in `DeConto2016_experiments.csv` but no longer shown in the hover.
- **DP4 → A:** climate_model is **CCSM4 (ocean); GENESIS v3 (atmosphere)**, and regional_climate_model is **RegCM3**.

**Data note:** member 1 has a LIG value of 3.53 m in Data Set S1, just under the 3.6 m lower target. The member is kept as Kopp et al. 2017 give it; a rounded comparison may explain it.

What each decision weighed:

- **DP1. Corrected vs uncorrected runs.**
  - **Option A:** switch to the "(BC)" rows, which are DP16 as published, with the +3 °C Amundsen–Bellingshausen ocean adjustment.
  - **Option B:** keep the uncorrected sensitivity rows.
  - **Option C:** bundle both, as separate experiments with an ocean-correction column. They'd be 58 member runs per scenario from one model configuration, and the uncorrected half would double-weight DP16.
  - Recommendation: **A**.
- **DP3. Ocean forcing and sensitivity for a non-ISMIP6 study.**
  - **Ocean forcing:** "Open" (the group's own melt parameterization, as for ISMIP6 open-framework runs).
  - **Ocean sensitivity:** a new value "Study-specific (see experiment)", with OCFAC shown per run, rather than forcing OCFAC onto the ISMIP6 percentile labels.
  - Recommendation: both as stated.
- **DP4. Climate model.**
  - **Option A:** a combined label, "CCSM4 (ocean); GENESIS v3 (atmosphere)", so DP16 doesn't pool with ISMIP6 runs whose atmosphere *and* ocean both come from CCSM4.
  - **Option B:** "CCSM4" alone, the only CMIP model involved.
  - Recommendation: **A**. DP16's retreat is driven largely by atmospheric meltwater (hydrofracturing), which isn't CCSM4's.
