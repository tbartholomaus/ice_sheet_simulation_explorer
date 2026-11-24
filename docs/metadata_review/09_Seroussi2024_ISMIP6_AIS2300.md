# Metadata review 9: Seroussi et al. (2024), ISMIP6 Antarctica 2300 (new study)

**Paper:** Seroussi, H., Pelle, T., Lipscomb, W. H., et al. (2024). Evolution of the Antarctic Ice Sheet over the next three centuries from an ISMIP6 model ensemble. *Earth's Future* 12, e2024EF004561. https://doi.org/10.1029/2024EF004561 (open access, CC BY 4.0). The paper was read in full, including Tables 1–4 and Appendix C. Citations to it below are "S24".

**Data:** "Scalars computed from two dimensional fields for this study ... are permanently available on Zenodo: https://doi.org/10.5281/zenodo.10528582" (CC BY 4.0). `ComputedScalars.zip` (127 MB) holds annual NetCDF scalars per experiment, variable and model.
- The app would use `ivaf`, "ice volume above floatation". Despite its name it is stored in **Gt** (units attribute "Gt"; values about 2.1 × 10⁷ Gt, Antarctica's ~58 m sea-level equivalent), so no density conversion is applied.
- It would cover 2015–2100 of each projection, joined to the model's historical run where one exists (32 of 43 configurations).

**Rules applied:** those of reviews 1–8 (`encoding_rules.csv`).

This is a **new study**, not yet in the app. The decisions below shape how it's added.

---

## 1. Summary

| # | Severity | Finding |
|---|---|---|
| 1 | **High** | **The archive's historical files have a mislabelled time axis.** Every `histAE` file starts at "2015" and counts forward: PIK_PISM's 165 values run "2015–2179". Each file's last value equals its projection's first value (2015) to within 10⁻⁷ Gt, and the counts fit each model's initial year in Table 3 (ILTS_SICOPOLIS, 1990: 25 values; LSCE_GRISLI, 1995: 21; ULB_Kori, 1950: 65). So the true years are 2015 − (n − 1) … 2015. The loader must rebuild them. |
| 2 | **High** | **Three experiments duplicate others through 2100.** The repeat-forcing runs repeat 2080–2100 forcing only after 2100 (S24 Sect. 2.5). expAE06 = expAE05, expAE08 = expAE03 and expAE09 = expAE04 through 2100, **exactly** (difference 0.0 Gt for every model). Keeping both would double-count those runs. Decision S3. |
| 3 | Medium | **S24 drops the standard/open split.** "Unlike Seroussi et al. (2020), we make no distinction between the ISMIP6 parameterizations and other approaches". Each group chose its melt scheme (Table 3). Decision S1. |
| 4 | Medium | **Ocean sensitivity is mostly unstated.** DOE_MALI gives it: "In our 4–20 km simulations, we used the median parameter value from the MeanAnt calibration, while for the 8–30 km simulations we ran separate ensembles with both the median and 95th percentile values." PIK_PISM uses Reese et al. (2023) coefficients. Other groups don't say. Decision S2. |
| 5 | Medium | **NORCE_CISM2's resolution conflicts between table and text.** Table 3 gives 8 km. Appendix C8 says "4 km (CISM2), 8 km (CISM3), 16 km (CISM4) or 32 km (CISM5)". Decision S4. |
| 6 | Medium | **IMAU_UFEMISM1–4 aren't distinguished anywhere.** Appendix C5 says "a medium-high resolution near the grounding line (either 30 km or 16 km depending on the model version)", but doesn't say which versions use which. Table 3 gives "30–200" for all four. Decision S4. |
| 7 | Low | **VUW_PISM's Table 3 code is SP, but its spin-up resets thickness.** Appendix C16: "an iterative nudging procedure, in which incremental grid refinement steps are employed that also include resetting of ice thicknesses to initial values". By rule D1 that is **Hybrid (thickness-nudged)**, as Seroussi 2020's VUW row is already encoded. |
| 8 | Low | **LSCE_GRISLI's sliding text says "power law (Weertman, 1957)".** But "LSCE_GRISLI2 corresponds exactly to the version used for (Seroussi et al., 2020)", where the reviewed encoding is **Linear** (m = 1). Both GRISLI rows are encoded Linear. |
| 9 | Low | **What was verified.** The loader's 329 runs (expAE06/08/09 left out) match Table 4's counts for every experiment, 43 configurations in all. Every main submission's 2015–2300 loss under expAE03 matches Table 5 to within about 3%, the app's are slightly higher: DC_ISSM 1.01 vs 0.98 m, ILTS_SICOPOLIS 4.39 vs 4.25, NCAR_CISM2 3.44 vs 3.33, PIK_PISM 2.05 vs 1.98, ULB_Kori1 1.74 vs 1.68. That is consistent with the app's fixed 362.5 Gt per mm, against the paper's ocean-density conversion. Every run reaches 2100. |
| 10 | Low | **43 configurations** (Table 3), each independent under your rule. Several are sensitivity variants: NORCE (11), VUW (10), UFEMISM (4), MALI (3). So S24 would count as 43 configurations, against 16 for Seroussi 2020. |

---

## 2. Experiments (Table 1)

| Experiment | GCM (verbatim) | Encoded GCM | Scenario | Forcing after 2100 | Ice-shelf collapse | In the app |
|---|---|---|---|---|---|---|
| Historical | — | — | — | — | No | joined to each projection |
| ctrlAE | — | — | control | — | No | not loaded (no other study's control runs are) |
| expAE01 | NorESM1-M | NorESM1-M | RCP2.6 | Repeat | No | yes |
| expAE02 | CCSM4 | CCSM4 | RCP8.5 | To 2300 | No | yes |
| expAE03 | HadGEM2 | HadGEM2-ES | RCP8.5 | To 2300 | No | yes |
| expAE04 | CESM2 | CESM2-WACCM | SSP5-8.5 | To 2300 | No | yes |
| expAE05 | UKESM | UKESM1-0-LL | SSP5-8.5 | To 2300 | No | yes |
| expAE06 | UKESM | UKESM1-0-LL | SSP5-8.5 | Repeat | No | S3 (= expAE05 to 2100) |
| expAE07 | NorESM1-M | NorESM1-M | RCP8.5 | Repeat | No | yes |
| expAE08 | HadGEM2 | HadGEM2-ES | RCP8.5 | Repeat | No | S3 (= expAE03 to 2100) |
| expAE09 | CESM2 | CESM2-WACCM | SSP5-8.5 | Repeat | No | S3 (= expAE04 to 2100) |
| expAE10 | UKESM | UKESM1-0-LL | SSP1-2.6 | To 2300 | No | yes |
| expAE11–14 | CCSM4, HadGEM2, CESM2, UKESM | as above | RCP8.5 / SSP5-8.5 | To 2300 | **Yes** | yes, excluded from the bias analysis (as Seroussi 2020's collapse runs) |

**Climate model identifiers.** S24 Sect. 2.1: "two CMIP5 models, the Community Climate System Model (CCSM4) and the Hadley Center Global Environment Model (HadGEM2-ES), and two CMIP6 models, the Community Earth System Model (CESM2-WACCM) and the UK Earth System Model (UKESM1-0-LL)". CESM2-WACCM is not the CESM2 (CAM6) of ISMIP6 2100, so it keeps its own identifier.

**Regional climate model:** "None (AOGCM output used directly)", as for Seroussi 2020. GCM surface mass balance anomalies are added to each model's reference climatology (Sect. 2.2).

---

## 3. Model configurations

**Model families.** Stress balance, resolution, initialization and melt scheme are from Table 3; sliding laws are from Appendix C. Quotes, by family:
- **DC_ISSM, UCSD_ISSM:** "a Budd sliding law (Budd et al., 1979)". DC_ISSM uses "the effective pressure equal to the ice overburden pressure, assuming the subglacial water is connected to the ocean". → **Budd**.
- **DOE_MALI:** "a power-law sliding relationship (Weertman, 1957) ... We found an exponent of 1/5 to give the best trade-off ... so we use this value in all forward simulations". → **Power law (Weertman)**, m = 5.
- **IGE_ElmerIce:** "a linear friction law". → **Linear**.
- **UTAS_ElmerIce:** "A linear Weertman sliding law (Weertman, 1957)". → **Linear**.
- **ILTS_SICOPOLIS:** "a Weertman–Budd-type sliding law with sub-melt sliding (Sato & Greve, 2012) and subglacial hydrology". → **Budd**.
- **IMAU_UFEMISM:** "a Budd-type sliding law". → **Budd**.
- **VUB_AISMPALEO:** "a Weertman relation inversely proportional to the height above buoyancy". → **Budd**, cascade step 2.
- **LSCE_GRISLI:** see finding 8. → **Linear**.
- **NCAR_CISM:** "Basal friction combines power-law and Coulomb behavior, following either Schoof (2005) (submission 1) or Zoet and Iverson (2020) (submission 2)". → **Coulomb-limited**.
- **NORCE_CISM:** "similar to that of Zoet and Iverson (2020) in all experiments". → **Coulomb-limited**.
- **UCM_Yelmo:** "a regularized-Coulomb law scaled by effective pressure". → **Coulomb-limited**.
- **PIK_PISM:** "A generalized power law (Schoof & Hindmarsh, 2010) ... The basal friction coefficient depends on the effective pressure and till friction angle". This is PISM's Mohr–Coulomb yield stress. → **Coulomb-limited**.
- **VUW_PISM:** pseudo-plastic, q = 0.4 (s1: 0.3, s2: 0.5). → **Coulomb-limited**.
- **ULB_Kori:** "a Weertman sliding law, with sliding exponent m = 3". → **Power law (Weertman)**.
- **UNN_Úa:** "here Weertman sliding law was used". The exponent isn't stated. → **Power law (Weertman)**, medium confidence.

**Initialization notes:**
- **PIK_PISM:** a 25,000-year spin-up ensemble. "The ensemble member performing best in the aggregated score was chosen as the initial state". That is selection, not nudging, so it's **Spin-up / equilibrium** (Table 3: SP).
- **IGE_ElmerIce:** the inverted friction "is reduced by 10% all over Antarctica" to "match the rate of ice-sheet mass change estimated by the IMBIE team". It is still velocity inversion, so **Data assimilation** (DA). The IMBIE tuning is noted.

**All 43 configurations,** shown as source code → encoding. Ocean forcing and sensitivity follow the recommended options of S1 and S2.

| Configuration | Stress balance | Resolution (km) | Initialization | Initial year | Sliding law | Melt scheme → ocean forcing | Ocean sensitivity |
|---|---|---|---|---|---|---|---|
| DC_ISSM | SSA → Shallow | 2–50 → ≤5 km | DA → Data assimilation (velocity inversion) | 2015 | Budd | Quad. Non-local → Standard | Medium |
| DOE_MALI_4km | HO → Higher-order | 4–20 → ≤5 km | DA+ → Data assimilation (velocity inversion) | 2000 | Power law (Weertman) | Quad. Non-local → Standard | Medium |
| DOE_MALI_8km_Ant95 | HO → Higher-order | 8–30 → >5–10 km | DA+ → Data assimilation (velocity inversion) | 2000 | Power law (Weertman) | Quad. Non-local → Standard | High |
| DOE_MALI_8km_AntMean | HO → Higher-order | 8–30 → >5–10 km | DA+ → Data assimilation (velocity inversion) | 2000 | Power law (Weertman) | Quad. Non-local → Standard | Medium |
| IGE_ElmerIce | SSA → Shallow | 1–50 → ≤5 km | DA → Data assimilation (velocity inversion) | 1995 | Linear | PICO → Open | Medium |
| ILTS_SICOPOLIS | Hybrid → Shallow | 8 → >5–10 km | SP+ → Hybrid (thickness-nudged) | 1990 | Budd | Quad. Non-local → Standard | Medium |
| IMAU_UFEMISM1 | Hybrid → Shallow | 30–200 → >20 km | SP+ → Hybrid (thickness-nudged) | 2014 | Budd | Quad. Local → Standard | Medium |
| IMAU_UFEMISM2 | Hybrid → Shallow | 30–200 → >20 km | SP+ → Hybrid (thickness-nudged) | 2014 | Budd | Quad. Local → Standard | Medium |
| IMAU_UFEMISM3 | Hybrid → Shallow | 30–200 → >20 km | SP+ → Hybrid (thickness-nudged) | 2014 | Budd | Quad. Local → Standard | Medium |
| IMAU_UFEMISM4 | Hybrid → Shallow | 30–200 → >20 km | SP+ → Hybrid (thickness-nudged) | 2014 | Budd | Quad. Local → Standard | Medium |
| LSCE_GRISLI | Hybrid → Shallow | 16 → >10–20 km | SP+ → Hybrid (thickness-nudged) | 1995 | Linear | Quad. Non-local → Standard | Medium |
| LSCE_GRISLI2 | Hybrid → Shallow | 16 → >10–20 km | SP+ → Hybrid (thickness-nudged) | 1995 | Linear | Quad. Non-local → Standard | Medium |
| NCAR_CISM1 | L1L2 → Higher-order | 4 → ≤5 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Non-local → Standard | Medium |
| NCAR_CISM2 | L1L2 → Higher-order | 4 → ≤5 km | SP+ → Hybrid (thickness-nudged) | 2015 | Coulomb-limited | Quad. Local → Standard | Medium |
| NORCE_CISM2 | L1L2 → Higher-order | 8 → >5–10 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Non-local Slope → Standard | Medium |
| NORCE_CISM3 | L1L2 → Higher-order | 8 → >5–10 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Non-local Slope → Standard | Medium |
| NORCE_CISM3_nonlocal | L1L2 → Higher-order | 8 → >5–10 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Non-local → Standard | Medium |
| NORCE_CISM3_local | L1L2 → Higher-order | 8 → >5–10 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Local → Standard | Medium |
| NORCE_CISM4 | L1L2 → Higher-order | 16 → >10–20 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Non-local Slope → Standard | Medium |
| NORCE_CISM4_nonlocal | L1L2 → Higher-order | 16 → >10–20 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Non-local → Standard | Medium |
| NORCE_CISM4_local | L1L2 → Higher-order | 16 → >10–20 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Local → Standard | Medium |
| NORCE_CISM4_JRA | L1L2 → Higher-order | 16 → >10–20 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Non-local Slope → Standard | Medium |
| NORCE_CISM5 | L1L2 → Higher-order | 32 → >20 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Non-local Slope → Standard | Medium |
| NORCE_CISM5_nonlocal | L1L2 → Higher-order | 32 → >20 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Non-local → Standard | Medium |
| NORCE_CISM5_local | L1L2 → Higher-order | 32 → >20 km | SP+ → Hybrid (thickness-nudged) | 1995 | Coulomb-limited | Quad. Local → Standard | Medium |
| PIK_PISM | Hybrid → Shallow | 8 → >5–10 km | SP → Spin-up / equilibrium | 1850 | Coulomb-limited | PICO → Open | Medium |
| UCM_Yelmo | L1L2 → Higher-order | 16 → >10–20 km | SP+ → Hybrid (thickness-nudged) | 1990 | Coulomb-limited | Quad. Non-local → Standard | Medium |
| UCSD_ISSM | SSA → Shallow | 3–50 → ≤5 km | DA → Data assimilation (velocity inversion) | 2007 | Budd | PICOP → Open | Medium |
| ULB_Kori1 | Hybrid → Shallow | 16 → >10–20 km | DA* → Hybrid (thickness-nudged) | 1950 | Power law (Weertman) | PICO → Open | Medium |
| ULB_Kori2 | Hybrid → Shallow | 16 → >10–20 km | DA* → Hybrid (thickness-nudged) | 1950 | Power law (Weertman) | Quad. Non-local → Standard | Medium |
| UNN_Úa | SSA → Shallow | 1–40 → ≤5 km | DA → Data assimilation (velocity inversion) | 2000 | Power law (Weertman) | Quad. Local → Standard | Medium |
| UTAS_ElmerIce | SSA → Shallow | 1–25 → ≤5 km | DA → Data assimilation (velocity inversion) | 1995 | Linear | Quad. Local → Standard | Medium |
| VUB_AISMPALEO | SIA + SSA → Shallow | 20 → >10–20 km | SP → Spin-up / equilibrium | 2000 | Budd | Quad. Non-local → Standard | Medium |
| VUW_PISM1 | Hybrid → Shallow | 8 → >5–10 km | SP → Hybrid (thickness-nudged) | 2015 | Coulomb-limited | Lin → Open | Medium |
| VUW_PISM1_s1 | Hybrid → Shallow | 8 → >5–10 km | SP → Hybrid (thickness-nudged) | 2015 | Coulomb-limited | Lin → Open | Medium |
| VUW_PISM1_s2 | Hybrid → Shallow | 8 → >5–10 km | SP → Hybrid (thickness-nudged) | 2015 | Coulomb-limited | Lin → Open | Medium |
| VUW_PISM1_s3 | Hybrid → Shallow | 8 → >5–10 km | SP → Hybrid (thickness-nudged) | 2015 | Coulomb-limited | Lin → Open | Medium |
| VUW_PISM1_s4 | Hybrid → Shallow | 8 → >5–10 km | SP → Hybrid (thickness-nudged) | 2015 | Coulomb-limited | Lin → Open | Medium |
| VUW_PISM2 | Hybrid → Shallow | 8 → >5–10 km | SP → Hybrid (thickness-nudged) | 2015 | Coulomb-limited | Lin → Open | Medium |
| VUW_PISM2_s1 | Hybrid → Shallow | 8 → >5–10 km | SP → Hybrid (thickness-nudged) | 2015 | Coulomb-limited | Lin → Open | Medium |
| VUW_PISM2_s2 | Hybrid → Shallow | 8 → >5–10 km | SP → Hybrid (thickness-nudged) | 2015 | Coulomb-limited | Lin → Open | Medium |
| VUW_PISM2_s3 | Hybrid → Shallow | 8 → >5–10 km | SP → Hybrid (thickness-nudged) | 2015 | Coulomb-limited | Lin → Open | Medium |
| VUW_PISM2_s4 | Hybrid → Shallow | 8 → >5–10 km | SP → Hybrid (thickness-nudged) | 2015 | Coulomb-limited | Lin → Open | Medium |
**Archive names.** The archive names some configurations differently from Table 3. These were matched to Table 3 by resolution, melt scheme and count:
- `NORCE_CISM4-MAR364-ERA-t1-nonlocal` = NORCE_CISM4_nonlocal
- `NORCE_CISM4-MAR364-JRA-t1` = NORCE_CISM4_JRA
- `ULB_fETISh-KoriBU1` = ULB_Kori1
- `UNN_Ua` = UNN_Úa

**Groups.** Group names get the prefix "Seroussi2024-" (e.g. `Seroussi2024-NCAR`), as Goelzer 2025's do. That keeps these configurations apart from Seroussi 2020's, even where a group and model name repeat (VUB AISMPALEO).

---

## 4. Decisions (applied 2026-11-24)

The user chose the recommended option of every decision.
- **S1:** ocean forcing by melt scheme, encoded per configuration (`ocean_forcing_encoded`). Standard for the ISMIP6 quadratic schemes (29), Open for PICO, PICOP and linear (14).
- **S2:** ocean sensitivity High for DOE_MALI_8km_Ant95 and Medium for the rest, encoded per configuration (`ocean_sensitivity_encoded`). `metadata.py` now supports a per-configuration sensitivity, as it already did for ocean forcing.
- **S3:** expAE06, 08 and 09 are left out (the loader's default). The app loads 329 runs.
- **S4:** NORCE_CISM2 at **4 km** (≤5 km), following Appendix C8. IMAU_UFEMISM1–4 at Table 3's 30 km, with a note.
- **S5:** every configuration counts on its own. In the app, that means S24 is weighted as an **ISMIP6 study** (`data.ISMIP6_STUDIES`): each of its runs counts once, like Seroussi 2020's. Other studies are each scaled to one typical ISMIP6 configuration. The notebook adds S24 as an extra source, which weighs it as one study, unlike the app.

**Also applied:**
- **Metadata files:** `Seroussi2024_models.csv` (43 rows), `Seroussi2024_experiments.csv` (Table 1 verbatim), and their column map.
- **Data and app:** the bundled data are `external_sources_seroussi2024_ais.csv.gz`. In the app the study is labelled "Seroussi 2024 (ISMIP6 AIS 2300)", listed in publication order after the two 2020 ISMIP6 ensembles (between Coulon 2024 and Rahlves 2025), coloured #d62728 and unchecked by default.

What each decision weighed:

- **S1. Ocean forcing.**
  - **Option A:** by melt scheme. ISMIP6 quadratic local, non-local or non-local-slope (Jourdain et al. 2020) → **Standard** (29 configurations); PICO, PICOP or linear → **Open** (14). This applies the existing rule's definition of Standard, "the quadratic sub-shelf melt parameterization (Jourdain et al. 2020)".
  - **Option B:** **Open** for all 43, following S24's statement that the distinction no longer applies.
  - Recommendation: **A**. It keeps the field meaningful across studies.
- **S2. Ocean sensitivity.**
  - **Option A:** **High** for DOE_MALI_8km_Ant95 (MeanAnt 95th percentile); **Medium** for every other configuration, as "the group's own calibration" under the existing rule. MALI's 4 km and AntMean runs are the MeanAnt median.
  - **Option B:** Medium only for MALI's median runs; "Study-specific (see experiment)" for the rest, since their calibrations aren't reported.
  - Recommendation: **A**.
- **S3. Duplicated experiments.**
  - **Option A:** leave out expAE06, expAE08 and expAE09, which are identical to expAE05, 03 and 04 through 2100.
  - **Option B:** keep them, double-counting those runs.
  - Recommendation: **A**.
- **S4. Unresolved resolutions.**
  - **NORCE_CISM2:** **Option A**, follow Appendix C8 (4 km, ≤5 km); **Option B**, follow Table 3 (8 km, >5–10 km). Recommendation: **A**. The appendix's sequence 4/8/16/32 km for CISM2–5 is internally consistent, while Table 3 lists CISM2 and CISM3 identically.
  - **IMAU_UFEMISM1–4:** encode all four at Table 3's 30 km (>20 km), with a note that some versions may use 16 km. Or ask the authors which versions do. Recommendation: **30 km with the note.**
- **S5. Weight.** Every configuration counts on its own, under your earlier rule, so S24 weighs as 43 ISMIP6 configurations.
  - **Option A:** keep that.
  - **Option B:** count each group's main submission (Table 3 asterisk) fully and its variants at reduced weight.
  - Recommendation: **A**, for consistency with your rule.
