# Metadata review 6: Coulon et al. (2024), Kori-ULB Antarctic ensemble

**Paper:** Coulon, V., Klose, A. K., Kittel, C., Edwards, T., Turner, F., Winkelmann, R. and Pattyn, F. (2024). Disentangling the drivers of future Antarctic ice loss with a historically calibrated ice-sheet model. *The Cryosphere* 18, 653–681. https://doi.org/10.5194/tc-18-653-2024 (open access). The main text, Table 1, Table 2 and Appendices A–C were read. Citations to the paper below are "C24".

**Data:** Coulon et al. (2023), Zenodo, https://doi.org/10.5281/zenodo.8398771 (version record 10812218). Its README.txt was read. The files used are:
- `HIST_ENSEMBLE_DATA.mat`: 1950–2014.
- `SSP585_ENSEMBLE_DATA.mat`: 2015 onward.
- `LHSensemble.mat`: the 100 × 9 parameter table.

The app uses all 100 members, 1950–2100, under SSP5-8.5.

**Rules applied:** those of reviews 1–5 (`encoding_rules.csv`).

---

## 1. Summary

| # | Severity | Finding |
|---|---|---|
| 1 | **High** | **The initialization is mis-encoded** as "Data assimilation (unreviewed)". C24 nudges thickness, not velocities (quote below the table). That is **Hybrid (thickness-nudged)** by rule D1. Seroussi 2020's ULB f.ETISh rows (code DA\*) are already Hybrid. |
| 2 | **High** | The other model metadata was entered without a source. C24 Appendix A gives hybrid SIA/SSA (**Shallow**), a Weertman law with **m = 3** (**Power law (Weertman)**) and **16 km** (**>10–20 km**). |
| 3 | **High** | **The loader drops one year of mass change at 2014/2015.** Each archive file starts at zero in its first year. The values are the state at the start of each year: SLC(t+1) − SLC(t) matches the archive's annual net mass balance for year t (r = 0.9999). The loader sets SSP 2015 equal to the historical 2014 value, so each member's 2014 change is lost (median 0.36 mm SLE, 130 Gt). Decision C4. |
| 4 | Medium | **The ice model was labelled differently in the two studies that run this code**: "Kori-ULB/f.ETISh" here and "f.ETISh (Kori-ULB)" in Seroussi 2020, so grouping by ice sheet model split one code base into two. C24: "Kori-ULB (previously f.ETISh)". Decision C1. |
| 5 | Medium | Ocean forcing and sensitivity aren't encoded. Each member draws one of five melt schemes and a continuous heat-exchange coefficient for it. Decision C2. |
| 6 | Medium | `protocol` reads "Not applicable". The members vary nine parameters, the 1950–2014 hindcast is forced by NorESM1-M (CMIP5), and the paper calibrates the ensemble against IMBIE. None of this is shown (§2). |
| 7 | Low | The regional climate model isn't encoded, and each member's present-day climatology (MARv3.11 or RACMO2.3p2) isn't shown (§2). |
| 8 | Low | The app uses the ensemble unweighted, while C24's projections weight members by their fit to IMBIE 1992–2011. Decision C5. |

**What finding 1 rests on.** C24 Appendix A: "Ice-sheet initial conditions and basal sliding coefficients are provided by an inverse simulation following Pollard and DeConto (2012b) ... basal sliding coefficients under grounded ice and sub-shelf melt rates under floating ice (Bernales et al., 2017) are adjusted iteratively to reduce the misfit with observed ice thickness".

**Effect of finding 3 on window rates** (member medians):

| Window | Change in rate |
|---|---|
| 2015–2023 (default) | none: only SSP5-8.5 values are fitted |
| 2010–2020 | −18 Gt/yr (−152 → −175) |
| 2012–2017 | −34 Gt/yr |
| 2005–2015 | −6 Gt/yr |

**What was verified.**
- **Parameter table.** The bundled data's per-member GCM and melt scheme match `LHSensemble.mat` columns 9 and 7 for all 100 members.
  - The four GCMs and five melt schemes each appear in 25 and 20 members.
  - The two atmospheric climatologies and the two ocean climatologies each appear in 50 members.
  - The coefficient ranges match Table 1. For example, γ0 for the ISMIP6 non-local scheme spans 1.07–3.86 × 10⁴ m/yr against Table 1's 1–4 × 10⁴.
- **Calibration periods.** The IMBIE calibration data (Table 2) cover 1992–2011 ("we select IMBIE estimates that correspond to our historical simulation, specifically those prior to 2015"). They don't overlap the app's default 2015–2023 window.

---

## 2. Experiments

- **Scenario:** SSP5-8.5 from 2015. The archive also has SSP1-2.6 and a constant-climate control, which the app doesn't use. Correct.
- **Climate forcing.**
  - "Hindcasts of the behaviour of the AIS over the period 1950–2014 are produced using changes in oceanic and atmospheric boundary conditions derived from the CMIP5 climate model NorESM1-M".
  - "As of the year 2015, changes in atmospheric and oceanic properties derived from a subset of CMIP6 climate models (MRI-ESM2-0, IPSL-CM6A-LR, CESM2-WACCM, and UKESM1-0-LL) are used".
  - Encoding: `climate_model` is the member's CMIP6 model, which forces the 2015–2023 window; the identifiers are already official CMIP6 IDs. The hover will add "1950–2014 hindcast: NorESM1-M".
- **Regional climate model.**
  - C24: "we compute future Antarctic SMB based on GCM outputs directly (similar to ... Nowicki et al., 2020; Seroussi et al., 2020)". GCM anomalies relative to 1995–2014 are added to a present-day climatology from MARv3.11 or RACMO2.3p2.
  - This is the ISMIP6-AIS approach, which is encoded "None (AOGCM output used directly)". Encoding: the same, with the member's climatology (MARv3.11 or RACMO2.3p2) shown in the hover.
- **Ocean forcing.** Each member uses one of five melt schemes:
  - PICO
  - Plume
  - Local quadratic
  - ISMIP6 non-local quadratic
  - ISMIP6 non-local quadratic with slope dependency

  Each member also has its own coefficient for its scheme, and uses either the Jourdain et al. (2020) or the Schmidtko et al. (2014) ocean climatology.

  For the two ISMIP6 schemes, γ0 is sampled over 1–4 × 10⁴ (or 1–4 × 10⁶) m/yr, with "the sector temperature corrections derived for the MeanAnt median γ0 values" applied at every γ0, and also with the Schmidtko climatology. So even these members don't follow the ISMIP6 standard protocol, which fixes γ0 at the MeanAnt 5th, 50th or 95th percentile. Decision C2.
- **Ice-shelf collapse:** none imposed (the hydrofracturing runs are a separate archive file, not used). `ice_shelf_fracture` = No.

---

## 3. Model configuration (Kori-ULB)

| Characteristic | Sources | Encoding |
|---|---|---|
| Ice model | "the Kori-ULB model", "Kori-ULB (previously f.ETISh)" (C24 Appendix A); README: "Kori-ULB ice-sheet model v0.9" | Decision C1 |
| Stress balance | "a vertically integrated, thermomechanical, hybrid ice-sheet–ice-shelf model ... a combination of the shallow-ice (SIA) and shallow-shelf (SSA) approximations for grounded ice, while only the shallow-shelf approximation is applied for floating ice shelves" | **Shallow** |
| Sliding law | "Basal sliding is introduced as a Weertman sliding law, i.e. vb = −Ab \|τb\|^(m−1) τb ... Ab the basal sliding coefficient, whose values are inferred following the nudging method of Pollard and DeConto (2012b); and m = 3 a sliding exponent". No effective pressure. | **Power law (Weertman)**, m = 3 (cascade step 4) |
| Initialization | Nudged inverse simulation (quoted under finding 1). "different initial states are only produced for each atmospheric present-day climatology"; then "a short 10-year relaxation". "Our ice-sheet simulations begin in the year 1950". | **Hybrid (thickness-nudged)** (D1), initial year 1950 |
| Resolution | "All simulations are performed at a spatial resolution of 16 km" | **>10–20 km** |

---

## 4. Recommended encoding

**Model row `Coulon2024 / Kori-ULB`** (one row, decision C3):

| Field | Value |
|---|---|
| ice_model | Kori-ULB/f.ETISh (C1) |
| stress_balance | Shallow |
| resolution | >10–20 km (16 km) |
| initialization | Hybrid (thickness-nudged), 1950 |
| sliding_law | Power law (Weertman), m = 3, no effective pressure, high confidence |

**Experiments:**

| Field | Value |
|---|---|
| climate_model | the member's CMIP6 GCM (unchanged) |
| scenario | SSP5-8.5 |
| regional_climate_model | None (AOGCM output used directly) |
| ocean_forcing | Open (C2) |
| ocean_sensitivity | Study-specific (see experiment) (C2) |
| ice_shelf_fracture | No |
| hover | melt scheme and coefficient, atmospheric and ocean climatologies, "1950–2014 hindcast: NorESM1-M", and the calibration note (C5) |

---

## 5. Decisions (applied 2026-11-23)

The user chose the recommended option for C2–C5. For C1 the user chose one label for both studies, **"Kori-ULB/f.ETISh"**.
- **C1:** both Coulon 2024 and Seroussi 2020's ULB rows (fETISh_16km, fETISh_32km) are encoded **Kori-ULB/f.ETISh**.
- **C2:** ocean forcing is **Open** and ocean sensitivity **Study-specific (see experiment)** for all 100 members. The hover shows each member's melt scheme and coefficient.
- **C3:** one model row.
- **C4:** the loader (`utilities/external_sources.py`) adds each member's 2014 change when joining the files, and `external_sources_coulon2024_ais.csv.gz` was regenerated. Values through 2014 are unchanged; from 2015 each member shifts by its own 2014 change (median −130 Gt).
- **C5:** the ensemble stays unweighted, and the hover says the paper weights members by fit to IMBIE 1992–2011.

`Coulon2024_experiments.csv` now transcribes all nine sampled parameters from `LHSensemble.mat`. The surface-melt parameters are kept in the file but not shown in the hover.

What each decision weighed:

- **C1. Ice model label.**
  - **Option A:** label Coulon's runs "f.ETISh (Kori-ULB)", the same as Seroussi 2020's ULB rows, so grouping by ice sheet model keeps the code base together, as rule G1 does for Elmer/Ice and CISM.
  - **Option B:** relabel both studies "Kori-ULB (f.ETISh)", the model's current name.
  - **Option C:** keep them separate.
  - Recommendation: **A**.
- **C2. Ocean forcing and sensitivity.**
  - **Option A:** ocean forcing **Open** for all 100 members (the group's own choice of schemes, with sampled coefficients and its own ocean climatologies). Ocean sensitivity **"Study-specific (see experiment)"**, with the scheme and coefficient in the hover.
  - **Option B:** the 40 ISMIP6 non-local members get ocean forcing "Standard", with their γ0 binned against the MeanAnt percentiles.
  - Recommendation: **A**, because none of the members follows the ISMIP6 standard protocol (§2).
- **C3. One configuration or two.**
  - The 100 members share one model version but start from two initial states, one per atmospheric climatology (MAR or RACMO). The model-level classes are the same for both.
  - **Option A:** one row, as for the DeConto parametric ensembles. Coulon then weighs as one ISMIP6 configuration.
  - **Option B:** two rows ("Kori-ULB (MAR)" and "Kori-ULB (RACMO)"), which doubles Coulon's weight.
  - **Correction (2026-11-24):** the weight statements above were wrong. In the app (`analysis.run_weights`), a non-ISMIP6 study as a whole weighs as one typical ISMIP6 configuration, however many configuration rows it has. Either option leaves Coulon's weight unchanged; the choice affects only how its runs group by model metadata.
  - Recommendation: **A**.
- **C4. Fix the 2014 gap.**
  - Restore each member's 2014 change when joining the files. Use the archive's net mass balance for 2014 (`mbcomp_ensemble` component 13), scaled by the member's own ratio of SLC step to net mass balance over 2004–2013 (median 1.03).
  - The default 2015–2023 rate is unchanged; windows spanning 2014 become correct.
  - Recommendation: **yes**.
- **C5. Weighting.**
  - **Option A:** keep the unweighted 100-member (prior) ensemble, with the hover noting that C24's projections weight members by fit to IMBIE 1992–2011.
  - **Option B:** apply C24's Bayesian weights within the study.
  - Recommendation: **A**. It's the ensemble as run, and the calibration doesn't touch the 2015–2023 window.
