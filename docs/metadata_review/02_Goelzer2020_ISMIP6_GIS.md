# Metadata review 2: Goelzer et al. (2020), ISMIP6-Greenland

**Paper:** Goelzer, H., et al. (2020). The future sea-level contribution of the Greenland ice sheet: a multi-model ensemble study of ISMIP6. *The Cryosphere* 14, 3071–3096. https://doi.org/10.5194/tc-14-3071-2020 (CC BY 4.0).

**Data:** `dash_app/data/ismip6_gis_ctrl.csv.gz`: 165 runs, 21 model configurations, built from the archive `v7_CMIP5_pub` ("v7_pub – July 15 2020 – same as v7 but renamed models to match publication").

**Rules applied:** the decisions from review 1 (tables first, appendix resolves ambiguity, conflicts flagged; three initialization classes; four-class sliding cascade; Shallow / Higher-order stress balance; binned finest resolution; every configuration independent). Where those rules don't settle a case, it is listed in §6 as a decision.

---

## 1. Summary

| # | Severity | Finding |
|---|---|---|
| 1 | **High** | **Four configurations ran the "standard" experiments in the *open* framework**, without the retreat parameterization (Table 3 footnote 1): BGC-BISICLES (all 9 experiments), UAF-PISM2, UCIJPL-ISSM2 and VUW-PISM. Ocean forcing is therefore a property of model × experiment for Greenland, not of the experiment alone. BGC-BISICLES also used its "own strategy" for the high/low ocean forcing of exp09/exp10 (footnote 2). (24 runs: 9 + 7 + 4 + 4) |
| 2 | **High** | **Sliding laws differ from the current encodings for most configurations.** Model-description papers resolve all 21 (§3). In particular: every ISSM configuration (AWI ×3, GSFC, JPL ×2, UCIJPL ×2) uses an effective-pressure (Budd) law; UAF-PISM is pseudo-plastic (q = 0.6); MUN-GSM1 is Coulomb-plastic on soft beds. |
| 3 | **High** | **Table A1 and Appendix A disagree on the stress balance of JPL-ISSM** (Table: HYB; App. A6: L1L2 over the whole domain) **and JPL-ISSMPALEO** (Table: SSA; App. A7 and Cuzzone et al. 2019: higher-order Blatter–Pattyn). Under the agreed scheme, L1L2 and HO are Higher-order, so the conflict changes the class. BGC-BISICLES (Table: SSA) uses the vertically integrated Schoof–Hindmarsh (2010) balance, which also keeps vertical shear in the rheology. Decision D3. |
| 4 | Medium | Initialization codes differ from Seroussi 2020's (DAv, DAs, DAi, NDs, NDm, CYC). They map onto the three classes as in §3b. Two cases need a decision: VUW-PISM (D1) and JPL-ISSMPALEO (D2). |
| 5 | Medium | **GIS ocean sensitivity is a different quantity** from AIS ocean sensitivity: Low/Medium/High are the 25th/50th/75th percentiles of the retreat parameter κ (Sect. 2; Slater et al. 2019, 2020), versus the 5th/50th/95th percentiles of the melt parameter γ0 for AIS. Decision D6. |
| 6 | Medium | **GCM names** differ in spelling from Seroussi 2020 for the same CMIP5 models (Table 1: "NorESM-M" (Table 3: "NorESM"), "CSIRO-Mk3.6", "ACCESS1-3"; AIS: "NorESM1-M", "CSIRO-MK3"). Without a common name the "GCM" factor can't pool across ice sheets. Decision D5. |
| 7 | Low | **Data gap, not metadata:** Table 3 lists MUN-GSM2 expa01–expa03, but the archive (and so the app) has no such runs. Every other configuration's experiments match Table 3 exactly. |
| 8 | Low | Small Table A1 vs Appendix A differences: IMAU initial year (Table: 1990; App. A5: steady state assigned to 1959, historical run 1960–2014; the archived series starts in 1961); UAF-PISM resolution (Table: 0.9 km; App. A11: regridded to 1 km for ISMIP6). Neither changes a class. |

**Mapping checks that passed:**
- The archive's model names match Goelzer's ("renamed models to match publication").
- The 165 runs match Table 3 apart from #7.
- The RCP assignment in `ismip6_gis_to_csv` (only exp07 is RCP2.6) matches Table 1.
- The bundled values are already control-run-corrected ("cr" files), as the paper's numbers are.

---

## 2. Experiments vs Table 1 (and Table 3)

Table 1 (verbatim): exp05 MIROC5 / 8.5 / Medium; exp06 NorESM-M / 8.5 / Medium; exp07 MIROC5 / 2.6 / Medium; exp08 HadGEM2-ES / 8.5 / Medium; exp09 MIROC5 / 8.5 / High; exp10 MIROC5 / 8.5 / Low; expa01 IPSL-CM5A-MR / 8.5 / Medium; expa02 CSIRO-Mk3.6 / 8.5 / Medium; expa03 ACCESS1-3 / 8.5 / Medium. Table 3 groups exp05–exp10 as "Core experiments" and expa01–expa03 as "Extensions". Experiments exp01–exp04 (open framework) are not in the archive.

The current `Goelzer2020_experiments.csv` (unreviewed) agrees on scenario and sensitivity. It differs in GCM spelling ("NorESM", "CSIRO-Mk3.6", "ACCESS1.3") and has no record of which runs were open-format.

---

## 3. Model configurations

### 3a. Sources used (beyond Goelzer 2020 Table A1 and Appendix A)

| Key | Source | Used for |
|---|---|---|
| R20 | Rückamp, Goelzer & Humbert 2020, TC 14, 3309 (the AWI ISMIP6-Greenland paper) | AWI-ISSM1–3 friction law, stress balance, resolution |
| L15 | Lee, Cornford & Payne 2015, Ann. Glaciol. 56(70), 129 | BGC-BISICLES friction law and stress balance (App. A2: initialization "the same as in initMIP-Greenland") |
| N23 | Nias, Nowicki, Felikson & Loomis 2023, JGR Earth Surface (NASA manuscript copy) | GSFC-ISSM friction law; it says it initializes "following the method used for the GSFC-ISSM contribution to ISMIP6" |
| F23 | Felikson et al. 2023, TC | GSFC-ISSM inversion following Morlighem et al. 2010 |
| C19 | Cuzzone et al. 2019, TC 13, 879 (cited by App. A7) | JPL-ISSMPALEO friction law and stress balance |
| Sch18 | Schlegel et al. 2018, TC 12, 3511 (review 1) | JPL ISSM friction law (Antarctic configuration, same group and code) |
| LAR | Levermann et al. 2020, LARMIP-2 (review 1) | UCI ISSM friction law (Antarctic configuration, same group and code) |
| dB15 | de Boer et al. 2015, TC 9, 881 | IMAU-ICE (SIA-only version) sliding |
| A16 | Aschwanden, Fahnestock & Truffer 2016, Nat. Commun. 7, 10524 (cited by App. A11 as the source of the initial state) | UAF-PISM sliding law and parameters |
| Q18, LC19 | Quiquet et al. 2018; Le clec'h et al. 2019 (GRISLI on Greenland) | LSCE-GRISLI friction law |
| L19 | Lipscomb et al. 2019, GMD 12, 387 (CISM v2.1) | NCAR-CISM sliding options |
| T25 | Tarasov et al. 2025, GMD 18, 9565 (GSM 25G) | MUN-GSM drag-law context |
| F15 | Fürst, Goelzer & Huybrechts 2015, TC 9, 1039 | VUB-GISM sliding law and stress balance |
| S19, ABU | initMIP-Greenland (Goelzer et al. 2018) and ABUMIP (Sun et al. 2020) | VUW/ARC-PISM pseudo-plastic sliding |

### 3b. Initialization: Table A1 codes mapped to the three classes

Table A1 codes: DAv = data assimilation of velocity; DAs = of surface elevation; DAi = of ice thickness; SP = spin-up; CYC = transient glacial cycle(s); NDm = nudging to ice mask; NDs = nudging to surface elevation. Footnote 7: "CYC/SP used only for the ice temperature".

| Code(s) | Class | Configurations |
|---|---|---|
| DAv | Data assimilation (velocity inversion) | AWI-ISSM1/2/3, BGC-BISICLES, GSFC-ISSM, JPL-ISSM, UCIJPL-ISSM1/2 |
| DAs, DAi, NDs (thickness/surface targets) | Hybrid (thickness-nudged) | ILTSPIK-SICOPOLIS1/2 (CYC/NDs), LSCE-GRISLI (SP/DAs), NCAR-CISM (SP/DAi), UAF-PISM1/2 (CYC/NDs, flux correction), VUB-GISM (CYC/DAi) |
| SP, CYC with NDm (ice mask only) | Spin-up / equilibrium | IMAU-IMAUICE1/2 (SP/NDm), MUN-GSM1/2 (CYC/NDm) |
| SP with NDs (thickness resets), per D1 | Hybrid (thickness-nudged) | VUW-PISM (SP/NDs) |
| **Needs a decision** | | JPL-ISSMPALEO (SP/DAv; D2) |

### 3c. Per configuration

"S" = stress balance (Table A1 code → class); "Res" = Table A1 min–max km → bin; "Init" = Table A1 code, initial year → class; sliding gives the source wording → class.

| Config | S | Res | Init | Sliding law (source) | Sliding class | Conf. |
|---|---|---|---|---|---|---|
| AWI-ISSM1 | HO → Higher-order (Blatter–Pattyn, R20) | 1–7.5 → ≤5 | DAv 1990 → DA | "basal drag τb follows a linear viscous law", τb = −k²N vb, N = ρi g h + min(0, ρw g zb) (R20 Eq. 1) | Budd (linear) | high |
| AWI-ISSM2 | HO → Higher-order | 1–7.5 → ≤5 | DAv 1990 → DA | as ISSM1 (R20; G1000 = AWI-ISSM2) | Budd | high |
| AWI-ISSM3 | HO → Higher-order | 0.75–7.5 → ≤5 | DAv 1990 → DA | as ISSM1 (R20; G750 = AWI-ISSM3) | Budd | high |
| BGC-BISICLES | SSA (Table A1) vs Schoof–Hindmarsh 2010 vertically integrated (L15) → **D3** | 1.2–4.8 → ≤5 | DAv 2000 → DA | "linear viscous friction law with coefficient C" (L15); effective-pressure term not readable (equation is an image) | Linear | medium |
| GSFC-ISSM | SSA → Shallow (N23: SSA) | 0.5–25 → ≤5 | DAv 2007 (+50-yr relaxation) → DA | "Budd friction law", τb = α²N ub, N = g(ρi H + ρw zb) (N23 Eq. 2) | Budd (linear) | high |
| ILTSPIK-SICOPOLIS1 | SIA → Shallow | 5 → ≤5 | CYC/NDs 1990 → Hybrid (topography nudging, "implied SMB", sliding coefficients tuned to velocities in 20 regions) | "Weertman–Budd-type sliding law with sub-melt sliding" (App. A4) | Budd | high |
| ILTSPIK-SICOPOLIS2 | HYB → Shallow | 5 → ≤5 | as SICOPOLIS1 → Hybrid | as SICOPOLIS1 | Budd | high |
| IMAU-IMAUICE1 | SIA → Shallow ("run in SIA mode") | 16 → >10–20 | SP/NDm 1990 → Spin-up | Not stated in Goelzer 2020. dB15: the SIA-only IMAU-ICE "only uses Weertman sliding"; form and exponent not given | Power law (Weertman) | **medium-low; D7** |
| IMAU-IMAUICE2 | SIA → Shallow | 8 → >5–10 | SP/NDm 1990 → Spin-up | as IMAUICE1 | Power law (Weertman) | **medium-low; D7** |
| JPL-ISSM | HYB (Table A1) vs L1L2 (App. A6) → **D3** | 0.25–15 → ≤5 | DAv 1979 (+4.2 kyr relaxation, 1840–1979 historical spin-up) → DA | Not stated in Goelzer 2020 or Schlegel et al. 2016's main text; inversion "following Morlighem et al. (2010)"; the same group's ISSM law is Budd, τb = −α²N vb (Sch18; C19) | Budd (linear) | medium-high |
| JPL-ISSMPALEO | SSA (Table A1) vs Blatter–Pattyn HO (App. A7, C19) → **D3** | 3–30 → ≤5 | SP/DAv 1979 → **D2** | τb = −k²N vb, k from velocity inversion, scaled through time with basal temperature (C19 Eqs. 2–3) | Budd (linear) | high |
| LSCE-GRISLI | HYB → Shallow | 5 → ≤5 | SP/DAs 1995 → Hybrid (iterative drag inversion to observed thickness, Le clec'h et al. 2019) | "power-law basal friction (Weertman, 1957)" (App. A8); GRISLI: m = 1, "linear viscous till", inverted β held constant (Q18, LC19) | Linear | medium-high |
| MUN-GSM1 | HYB → Shallow | 5–14 (0.25° × 0.125°) → ≤5 | CYC/NDm 1980 → Spin-up (two glacial cycles; margins nudged to 2000 CE) | "Coulomb-plastic soft-bed rheology" since 1500 CE; "power law basal drag formulation for hard bed" (App. A9) | Coulomb-limited | high |
| MUN-GSM2 | HYB → Shallow | 5–14 → ≤5 | CYC/NDm 1980 → Spin-up | "linear Weertman-type basal drag law" on soft beds; power law on hard beds (App. A9) | **D4** | high (form) |
| NCAR-CISM | HO → Higher-order (depth-integrated higher order, Goldberg 2011) | 4 → ≤5 | SP/DAi 1990 → Hybrid (friction coefficients nudged toward observed thickness during a 30 kyr spin-up) | "a basal-sliding law based on Schoof (2005)"; "no dependence of basal sliding on basal temperature or water pressure" (App. A10) | Coulomb-limited | medium-high |
| UAF-PISM1 | HYB → Shallow | 0.9 → ≤5 | CYC/NDs 2008 → Hybrid (glacial-cycle spin-up, then flux correction keeping the surface near observations) | "pseudoplastic power law", q = 0.6, u0 = 100 m/yr, Mohr–Coulomb yield stress with till friction angle by bed elevation, till-water effective pressure (A16, the "calibrated" runs used as the ISMIP6 initial state) | Coulomb-limited | high |
| UAF-PISM2 | HYB → Shallow | 0.9 → ≤5 | as PISM1 → Hybrid | as PISM1; ran in open format (Table 3 fn. 1) | Coulomb-limited | high |
| UCIJPL-ISSM1 | HO → Higher-order | 0.5–30 → ≤5 | DAv 2007 (+50-yr relaxation) → DA | Not stated in Goelzer 2020; "data assimilation (Morlighem et al., 2010)"; the same group's ISSM inversions use a Budd law (LAR A8; N23/F23 follow the same method) | Budd (linear) | medium-high |
| UCIJPL-ISSM2 | HO → Higher-order | 0.2–20 → ≤5 | DAv 2007 → DA | as ISSM1; moving front; ran in open format | Budd (linear) | medium-high |
| VUB-GISM | HO → Higher-order ("SR HO"; F15: multilayer longitudinal-stress approximation, LMLa) | 5 → ≤5 | CYC/DAi 1990 → Hybrid (assimilation of observed thickness, Le clec'h et al. 2019) | Weertman sliding "extended with a multiplier", vb = SBL (AS/H) τb³ (F15 Eq. 1) | Budd (m = 3, ∝ 1/H) | medium-high |
| VUW-PISM | HYB → Shallow | 2 → ≤5 | SP/NDs 2000 → Hybrid (D1) | Not stated in Goelzer 2020 (refers to Golledge et al. 2019); initMIP-Greenland ARC-PISM: "pseudo-plastic exponent q"; ABUMIP ARC-PISM: "Coulomb q = 0.75" | Coulomb-limited | medium |

**Ocean forcing by configuration:** Standard for all, except **Open** for BGC-BISICLES, UAF-PISM2, UCIJPL-ISSM2 and VUW-PISM (Table 3 footnote 1).

---

## 4. Implementation notes

- Ocean forcing needs a **model-level override** for Greenland: a `ocean_forcing_encoded` column in `Goelzer2020_models.csv`, used where it is filled in, otherwise the experiment's value. `metadata.py` would validate it like the other encoded fields.
- `Goelzer2020_experiments.csv` becomes a verbatim transcription of Table 1 (`Exp ID`, `GCM`, `RCP`, `Ocean sensitivity`) plus Table 3's Core/Extensions grouping. Value maps translate "8.5" → "RCP8.5", and the GCM names per D5.

---

## 5. Places where Goelzer 2020 disagrees with itself (or with the archive)

| Item | Table A1 / Table 3 | Elsewhere |
|---|---|---|
| JPL-ISSM stress balance | HYB | App. A6: "L1L2 ... is used over the entire domain" |
| JPL-ISSMPALEO stress balance | SSA | App. A7: "higher-order ice flow approximation of Blatter (1995) and Pattyn (2003)"; C19 the same |
| BGC-BISICLES stress balance | SSA | L15: vertically integrated stress balance of Schoof and Hindmarsh (2010), with vertical shear in the rheology |
| IMAU initial year | 1990 | App. A5: steady state assigned to 1959; historical run 1960–2014 (archive starts 1961) |
| UAF-PISM resolution | 0.9 km | App. A11: regridded to 1 km for ISMIP6 |
| MUN-GSM2 experiments | Table 3: expa01–03 run | Not in the archive |

---

## 6. Decisions needed

- ~~**D1. VUW-PISM initialization.**~~ **DECIDED (2026-11-23):** any spin-up nudged to observed ice thickness or surface elevation is Hybrid (thickness-nudged). VUW-PISM is Hybrid on both ice sheets; the AIS row was changed to match.
- ~~**D2. JPL-ISSMPALEO initialization.**~~ **DECIDED:** Data assimilation (velocity inversion), by its velocity-inverted friction field.
- ~~**D3. Stress-balance conflicts.**~~ **DECIDED:** follow the appendix and model papers. JPL-ISSM (L1L2; also Schlegel et al. 2016), JPL-ISSMPALEO (Blatter–Pattyn; also Cuzzone et al. 2019) and BGC-BISICLES (Schoof–Hindmarsh / L1L2-type) are Higher-order.
- ~~**D4. MUN-GSM2 sliding.**~~ **DECIDED:** Linear (the soft-bed law).
- ~~**D5. GCM names.**~~ **DECIDED:** encode official CMIP model identifiers through value maps; the verbatim spellings stay in the experiments files. Applied to Seroussi2020 (CSIRO-MK3 → CSIRO-Mk3-6-0) and Goelzer2020 (NorESM-M → NorESM1-M, CSIRO-Mk3.6 → CSIRO-Mk3-6-0); other studies at their review.
- ~~**D6. Ocean sensitivity.**~~ **DECIDED:** shared Low/Medium/High labels; encoding_rules.csv gives each ice sheet's meaning (AIS γ0 5/50/95th, GIS κ 25/50/75th percentiles).
- ~~**D7. IMAU-ICE sliding.**~~ **DECIDED:** Power law (Weertman), medium-low confidence, on the strength of de Boer et al. (2015).


---

## 7. Implementation (2026-11-23)

On branch `metadata-review-ais`, commit `1fc791b`:

- `dash_app/metadata/Goelzer2020_models.csv`: all 21 configurations reviewed, with source wording, source, encoded value and notes for each characteristic, plus a per-configuration `ocean_forcing_encoded` (Open for BGC-BISICLES, UAF-PISM2, UCIJPL-ISSM2 and VUW-PISM).
- `dash_app/metadata/Goelzer2020_experiments.csv`: Table 1 verbatim (with Table 3's Core/Extensions grouping). `experiment_columns.csv` maps "8.5" → "RCP8.5", "NorESM-M" → "NorESM1-M" and "CSIRO-Mk3.6" → "CSIRO-Mk3-6-0".
- `encoding_rules.csv`: the initialization rules now cover Goelzer's codes (D1). The ocean-forcing and ocean-sensitivity meanings cover both ice sheets (D6). A new rule row records the CMIP naming convention (D5).
- `metadata.py`: an optional model-level ocean-forcing override (validated), and a generic experiments key column.
- `Seroussi2020`: VUW-PISM → Hybrid (D1); CSIRO-MK3 → CSIRO-Mk3-6-0 (D5).
- Tests: 61 pass. The run table changed in exactly the expected rows: 165 GIS runs, the 4 AIS VUW runs (initialization) and the 14 AIS CSIRO runs (GCM name). No external study changed.
