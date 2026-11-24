# Metadata review 1: Seroussi et al. (2020), ISMIP6-Antarctica

**Paper:** Seroussi, H., et al. (2020). ISMIP6 Antarctica: a multi-model ensemble of the Antarctic ice sheet evolution over the 21st century. *The Cryosphere* 14, 3033–3070. https://doi.org/10.5194/tc-14-3033-2020 (CC BY 4.0). Read from a local PDF copy.

**Code reviewed:**
- `dash_app/data.py`: `ais_exp_meta`, `ism_meta`, `get_exp_meta`, `get_ism_meta`, `_build_runs`
- `analyze_slr_predictions_interactive.ipynb`, cell 10. It holds an identical copy of all three tables. I checked this programmatically, so every finding here applies to both.
- `utilities/data_loader.py`: `ismip6_ais_to_csv`
- `utilities/helper.py`: `ais_exp_dict`
- `dash_app/data/ismip6_ais.csv.gz`: 188 runs
- The source archive `ComputedScalarsPaper/` at the repository root (Zenodo 3940766)

**Ranking of evidence (agreed):** the paper's tables come first. The Appendix C model descriptions resolve ambiguity. Every disagreement between them is flagged below with its location.

## Decisions log

| Date | Decision | Applies to |
|---|---|---|
| 2026-11-23 | **Initialization: Option C, three classes by method** (§6.1). The verbatim paper code is stored alongside the class. | All papers |
| 2026-11-23 | **Sliding law: four classes assigned by a cascade** (§6.2). Take the first rule that applies: (1) basal stress capped by a yield/Coulomb limit → **Coulomb-limited**; (2) else friction depends on effective pressure or height above flotation → **Budd**; (3) else exponent = 1 → **Linear**; (4) else → **Power law (Weertman)**. Also store the exponent and an effective-pressure flag where stated. | All papers |
| 2026-11-23 | **Sliding laws are checked against the wider literature**: fill the "not stated" models from cited model-description papers, and cross-check every stated law against other papers by the same group/model. | All papers |
| 2026-11-23 | **Stress balance: two classes** (§6.3 Option C). **Shallow** = Hybrid, SIA+SSA, SSA; **Higher-order** = L1L2, HO (Blatter–Pattyn), Stokes. The verbatim paper code is stored alongside the class. | All papers |
| 2026-11-23 | **Resolution: binned finest resolution** (§6.4 Option B). Store `res_min_km` and `res_max_km`; group and analyze by `res_bin` from `res_min_km`: **≤5 km**, **>5–10 km**, **>10–20 km**, **>20 km**. Bin edges belong to the lower bin (e.g. 10 km → >5–10, 20 km → >10–20). | All papers |
| 2026-11-23 | **Imposed ice-shelf-collapse runs are excluded from the bias analysis** (§6.5). For ISMIP6-AIS that means exp11 and exp12 (14 runs; 174 remain). They are flagged by `ice_shelf_fracture = True` and left out of the rate-bias, low-bias-enrichment and ANOVA calculations. They stay in the time-series and 2100 views, labeled. When later papers have similar experiments, the same rule applies. | All papers |
| 2026-11-23 | **Ocean-sensitivity variants stay in the bias analysis** (§6.5). For ISMIP6-AIS: exp09 (High), exp10 (Low) and exp13 (PIGL), 37 runs. Each run's `ocean_sensitivity` value (Low/Medium/High/PIGL; "N/A (open)" for open runs) is stored on every run. | All papers |
| 2026-11-23 | **Ocean forcing and ocean sensitivity become two separate ANOVA/grouping factors** (supersedes "store only" above; also settles the melt-framework question in §6.5). `ocean_forcing` = Open / Standard (Table 1 "Ocean forcing"). `ocean_sensitivity` = Low / Medium / High / PIGL (Table 1 "Ocean sensitivity"). Open runs take **Medium**, as Table 1 lists them (tables-first rule). This also avoids an "N/A (open)" level that would just duplicate `ocean_forcing = Open` and make the ANOVA design rank-deficient. | All papers |
| 2026-11-23 | **Weighting stays per run** (§6.5): every ISMIP6 run counts 1, as `run_weights` does now. Model variants (ULB 16/32 km, IMAU 1/2, PIK 1/2) are *not* merged. In the AIS bias set (174 runs), ULB's two variants carry 38 runs (22 %) and AWI_PISM1 carries 19. | All papers |
| 2026-11-23 | **Every model configuration/version is independent.** Each distinct (Group, Model) is its own simulation configuration everywhere: its own metadata row, its own unit for weighting, and never merged with or borrowed from a sibling (e.g. PIK_PISM1 vs PIK_PISM2 are different model versions). Consequence for `run_weights`: the "typical ISMIP6 model" used to scale external studies should be runs per *configuration*, not runs per *institution*. For AIS that changes it from 188/13 groups = 14.5 to 188/16 configurations = 11.8 runs; for GIS from 165/13 = 12.7 to 165/21 = 7.9. Each external study then weighs as much as one typical ISMIP6 configuration. | All papers |
| 2026-11-23 | **Any spin-up nudged to observed ice thickness or surface elevation is "Hybrid (thickness-nudged)"**, whatever the mechanism: friction/melt/SMB adjustment, flux correction, or thickness resets (decided during review 2, D1). Consequence for AIS: **VUW_PISM moves from Spin-up / equilibrium to Hybrid** (its multistage spin-up resets ice thickness to observed values). Nudging to the ice mask only stays Spin-up. | All papers |
| 2026-11-23 | **Climate models are encoded with their official CMIP model identifiers** through value maps; verbatim spellings stay in the experiments files (review 2, D5). AIS: CSIRO-MK3 → CSIRO-Mk3-6-0. | All papers |
| 2026-11-23 | **Ocean sensitivity keeps shared Low/Medium/High labels across ice sheets**, with each ice sheet's meaning documented in encoding_rules.csv (review 2, D6). | All papers |

---

## 1. Summary

| # | Severity | Finding | Runs affected (of 188) |
|---|---|---|---|
| 1 | **High** | `expA1`–`expA8` are labeled **SSP** scenarios. They are **CMIP5 RCP** experiments (Table 1). `expA3`/`expA7` are also given the wrong forcing level: they are labeled SSP1-2.6 but are **RCP8.5**. | 55 (14 of them are RCP8.5 labeled as low-emission) |
| 2 | **High** | Initialization is **"Data assimilation" for IMAU_IMAUICE1, IMAU_IMAUICE2, LSCE_GRISLI and NCAR_CISM**. Table 3 codes them **Eq, SP, SP+ and SP+**. | 51 |
| 3 | **High** | Sliding law is **"Coulomb" for DOE_MALI** and **"Regularized Coulomb" for UTAS_ElmerIce**. Appendix C says both use a **linear** friction law. | 11 |
| 3b | **High** | *(Added after the literature check, §9.)* **IMAU_IMAUICE1/2 "Weertman" → plastic Mohr–Coulomb (q = 0)**; **UCIJPL_ISSM "Regularized Coulomb" → Budd**; **LSCE_GRISLI "Weertman" → effectively linear (m = 1)**. | 49 |
| 4 | **High** | `basal_melt_param` for exp09/exp10/exp13 is "PIGL medium / PIGL high / PIGL very high". Table 1 says **exp09 = High, exp10 = Low** (MeanAnt γ0, 95th/5th percentile) and **exp13 = PIGL**. | 37 |
| 5 | **High** | **Ice-shelf collapse (exp11, exp12) is not encoded at all.** These runs are pooled with every other run in all statistics. | 14 |
| 6 | Medium | Sliding law is "Weertman / Coulomb" for ULB fETISh. Appendix C10 describes a **power law only**. | 42 |
| 7 | Medium | *(Resolved in §9.)* **Sliding law is not stated in Seroussi 2020** for AWI_PISM, IMAU_IMAUICE1/2, JPL1_ISSM, UCIJPL_ISSM and VUW_PISM. The current values ("Pseudo-plastic", "Weertman", "Budd / Schoof", "Regularized Coulomb", "Pseudo-plastic (Coulomb)") come from somewhere else and have no citation. | 74 |
| 8 | Medium | **Open vs. standard melt framework** is only shown in hover text. It is not a grouping or ANOVA dimension. Seroussi et al. analyze the two separately because the framework "affects the results significantly" (Sect. 3). For open runs, `basal_melt_param` = "Standard" is wrong: each model uses its own scheme (Table 3). | 52 open runs |
| 9 | Medium | **Mapping:** `ULB/fETISh_16km` and `ULB/fETISh_32km` have no exact key. They resolve only through `get_ism_meta`'s substring fallback, which happens to return `("ULB","FETISH1")`. | 42 |
| 10 | Medium | **Structure:** `ism_meta` is keyed by (Group, Model) with no ice sheet. Five keys are **shared by AIS and GIS runs**: IMAU/IMAUICE1, IMAU/IMAUICE2, LSCE/GRISLI, NCAR/CISM, VUW/PISM. So fixing an AIS entry (e.g. #2) would silently relabel the GIS runs too. | 55 AIS runs, plus GIS runs |
| 11 | Low | `climate_model` "NorESM" should be **NorESM1-M** (Table 1). "CSIRO-MK3" matches the paper; its CMIP5 name is CSIRO-Mk3-6-0. | 74 |
| 12 | Low | The code comment cites "Seroussi et al. 2020 Table A1" for model characteristics. In the paper, Table A1 is the data request. Model characteristics are in **Table 3** and **Appendix C**. | – |
| 13 | Low | The data CSV has a correct `RCP` column, derived from the experiment number in `ismip6_ais_to_csv`. `_build_runs` ignores it and uses the incorrect lookup value, so scenario is defined in two places that disagree. | – |

**Mapping checks that passed:**
- Run inventory: the 188 (Group, Model, Exp) runs in `ismip6_ais.csv.gz` match Table 4's experiment matrix **exactly**, model by model and experiment by experiment. I checked against an image of the table because the PDF's rotated headers don't extract cleanly.
- The archive directory names match the CSV's Group/Model strings.
- `ais_exp_dict` (helper.py), which pairs each projection with `hist_open` or `hist_std`, matches Table 1's open/standard assignment for all 21 experiments.
- The archive netCDFs carry no model metadata in their attributes, so the lookup tables are the only source.

---

## 2. Experiment metadata (`ais_exp_meta`) vs. Table 1

Table 1 columns: AOGCM, scenario, ocean forcing (= melt framework), ocean sensitivity, ice-shelf fracture, tier.

| Exp | Current: climate_model / scenario / protocol / basal_melt_param | Table 1 | Verdict |
|---|---|---|---|
| exp01 | NorESM / RCP8.5 / Open / Standard | NorESM1-M, RCP8.5, Open, Medium, no fracture, Tier 1 | GCM name; melt param |
| exp02 | MIROC-ESM-CHEM / RCP8.5 / Open / Standard | MIROC-ESM-CHEM, RCP8.5, Open, Medium, no | melt param |
| exp03 | NorESM / RCP2.6 / Open / Standard | NorESM1-M, RCP2.6, Open, Medium, no | GCM name; melt param |
| exp04 | CCSM4 / RCP8.5 / Open / Standard | CCSM4, RCP8.5, Open, Medium, no | melt param |
| exp05 | NorESM / RCP8.5 / Standard / Standard | NorESM1-M, RCP8.5, Standard, Medium, no | GCM name |
| exp06 | MIROC-ESM-CHEM / RCP8.5 / Standard / Standard | MIROC-ESM-CHEM, RCP8.5, Standard, Medium, no | OK |
| exp07 | NorESM / RCP2.6 / Standard / Standard | NorESM1-M, RCP2.6, Standard, Medium, no | GCM name |
| exp08 | CCSM4 / RCP8.5 / Standard / Standard | CCSM4, RCP8.5, Standard, Medium, no | OK |
| exp09 | NorESM / RCP8.5 / Standard / **PIGL medium** | NorESM1-M, RCP8.5, Standard, **High** (MeanAnt 95th pct.), no | **Wrong sensitivity** |
| exp10 | NorESM / RCP8.5 / Standard / **PIGL high** | NorESM1-M, RCP8.5, Standard, **Low** (MeanAnt 5th pct.), no | **Wrong sensitivity** |
| exp11 | CCSM4 / RCP8.5 / Open / Standard | CCSM4, RCP8.5, Open, Medium, **fracture: Yes** | **Collapse not encoded**; melt param |
| exp12 | CCSM4 / RCP8.5 / Standard / Standard | CCSM4, RCP8.5, Standard, Medium, **fracture: Yes** | **Collapse not encoded** |
| exp13 | NorESM / RCP8.5 / Standard / **PIGL very high** | NorESM1-M, RCP8.5, Standard, **PIGL**, no | **Wrong sensitivity** |
| expA1 | HadGEM2-ES / **SSP5-8.5** / Open / Standard | HadGEM2-ES, **RCP8.5**, Open, Medium, no, Tier 2 | **Wrong scenario family**; melt param |
| expA2 | CSIRO-MK3 / **SSP5-8.5** / Open / Standard | CSIRO-MK3, **RCP8.5**, Open, Medium, no, Tier 2 | **Wrong scenario family**; melt param |
| expA3 | IPSL-CM5A-MR / **SSP1-2.6** / Open / Standard | IPSL-CM5A-MR, **RCP8.5**, Open, Medium, no, Tier 2 | **Wrong family and level** |
| expA4 | IPSL-CM5A-MR / **SSP1-2.6** / Open / Standard | IPSL-CM5A-MR, **RCP2.6**, Open, Medium, no, Tier 2 | **Wrong scenario family**; melt param |
| expA5 | HadGEM2-ES / **SSP5-8.5** / Standard / Standard | HadGEM2-ES, **RCP8.5**, Standard, Medium, no, Tier 2 | **Wrong scenario family** |
| expA6 | CSIRO-MK3 / **SSP5-8.5** / Standard / Standard | CSIRO-MK3, **RCP8.5**, Standard, Medium, no, Tier 2 | **Wrong scenario family** |
| expA7 | IPSL-CM5A-MR / **SSP1-2.6** / Standard / Standard | IPSL-CM5A-MR, **RCP8.5**, Standard, Medium, no, Tier 2 | **Wrong family and level** |
| expA8 | IPSL-CM5A-MR / **SSP1-2.6** / Standard / Standard | IPSL-CM5A-MR, **RCP2.6**, Standard, Medium, no, Tier 2 | **Wrong scenario family** |

**Notes**
- Table 1's caption says every projection is "based on CMIP5 AOGCMs". None of the ISMIP6-AIS runs is SSP-forced.
- In the app, the scenario mislabel also affects the "RCP/SSP composite" view: expA3/expA7 are pooled into "composite 2.6".
- "Ocean sensitivity" is defined only for the standard framework: the 5th/50th/95th percentiles of the MeanAnt γ0 distribution, or the PIGL calibration. For open runs, Table 1 lists "Medium", but the melt rate actually comes from the group's own parameterization. So the open-run value should be that scheme (from Table 3), not "Standard".

---

## 3. Model metadata (`ism_meta`) vs. Table 3 and Appendix C

### 3a. Existing fields

"Paper" gives the Table 3 initialization code and a paraphrase of the Appendix C description.

| Group_Model | Current ice_model / sliding_law / initialization | Paper: sliding law (App. C) | Paper: initialization (Table 3 code; App. C) | Verdict |
|---|---|---|---|---|
| AWI_PISM1 | PISM / Pseudo-plastic / Spin-up | **Not stated.** C1 says only that "the SSA solution acts as a sliding law". | **Eq**, 2005. Equilibrium-type spin-up under steady present-day climate: 200 kyr thermal, then 30 kyr with free geometry. | Sliding law unsupported by this paper. Init depends on the category decision. |
| DOE_MALI | MALI / **Coulomb** / Data assimilation | **Linear** friction law, spatially varying coefficient (C2) | **DA+**, 2015. Adjoint inversion of friction and shelf viscosity, then 99-yr relaxation. | **Sliding law wrong.** Init OK at 2-class level. |
| ILTS_PIK_SICOPOLIS | SICOPOLIS / Weertman / Spin-up | Weertman–Budd type with sub-melt sliding (C3) | **SP+**, 1990. 140 kyr paleo spin-up with topography nudged to present; sliding coefficient tuned per IMBIE basin to match velocities. | OK at the current level. Effective-pressure (Budd) dependence not recorded. |
| IMAU_IMAUICE1 | IMAU-ICE / Weertman / **Data assimilation** | **Not stated.** C4: SIA+SSA velocities added to represent sliding. | **Eq**, 1978. Thermodynamic, then dynamic steady state under constant present-day forcing. | **Init wrong.** Sliding law unsupported. |
| IMAU_IMAUICE2 | IMAU-ICE / Weertman / **Data assimilation** | **Not stated** (C4) | **SP**, 1978 in Table 3. Same steady state, then a forced run from 1900 (see §5). | **Init wrong.** Sliding law unsupported. |
| JPL1_ISSM | ISSM / **Budd / Schoof** / Data assimilation | **Not stated.** C5: friction coefficient inferred by inversion; law form not given. | **DA**, 2007. Inversion followed by a 2-yr relaxation (see §5). | Sliding law unsupported by this paper. Init OK. |
| LSCE_GRISLI | GRISLI / Weertman / **Data assimilation** | Power law (Weertman 1957) (C6) | **SP+**, 1995. Iterative inversion of the basal drag coefficient to match observed *thickness* with minimal drift. | **Init disagrees with Table 3.** Sliding OK. |
| NCAR_CISM | CISM / Regularized Coulomb / **Data assimilation** | Similar to Schoof (2005), combining power-law and Coulomb behavior (C7) | **SP+**, 1995. 30 kyr spin-up with friction and sub-shelf melt nudged toward observed thickness; the text calls it "a hybrid approach between assimilation and spin-up". | **Init disagrees with Table 3.** Sliding OK. |
| PIK_PISM1 | PISM / Pseudo-plastic / Spin-up | Power law with Mohr–Coulomb yield stress from till properties and effective pressure (C8) | **SP**, 1850. Described as an "equilibrium simulation". | OK |
| PIK_PISM2 | PISM / Pseudo-plastic / Spin-up | Same as PISM1 (C8) | **SP**, 2015 | OK |
| UCIJPL_ISSM | ISSM / **Regularized Coulomb** / Data assimilation | **Not stated.** C9: two-step inversion; law form not given. | **DA**, 2007 | Sliding law unsupported by this paper. Init OK. |
| ULB_fETISh_16km | Kori-ULB/f.ETISh / **Weertman / Coulomb** / Data assimilation | Power-law sliding, "power law coefficient of 2" (C10) | **DA\*** (geometry only), 2005. Iterative Pollard & DeConto (2012a) procedure fitting thickness and flow. | **"/ Coulomb" unsupported.** Init depends on the category decision. Key resolves only via fallback (#9). |
| ULB_fETISh_32km | (same, via fallback) | Same | **DA\***, 2005 | Same |
| UTAS_ElmerIce | Elmer/Ice / **Regularized Coulomb** / Data assimilation | **Linear** sliding relation (C11) | **DA**, 2015. Three inversions for the sliding coefficient, then a 10-yr surface relaxation (see §5). | **Sliding law wrong.** Init OK. |
| VUB_AISMPALEO | AISMPALEO / Weertman / Spin-up | Weertman, inversely proportional to height above buoyancy (C12) | **SP**, 2000. Free geometry run to steady state. | OK. Effective-pressure dependence not recorded. |
| VUW_PISM | PISM / Pseudo-plastic (Coulomb) / Spin-up | **Not stated.** C13 defers to Golledge et al. (2019). | **SP**, 2015 in Table 3. Multistage spin-up with thickness resets (nudging); historical 1950–2000 (see §5). | Sliding law unsupported by this paper. Init OK. |

The `ice_model` names are all fine. One small point: the paper calls ULB's model f.ETISh v1.3. The current label "Kori-ULB/f.ETISh" uses the later name, which is fine as long as it's used consistently with Coulon 2024.

### 3b. New fields: stress balance and resolution (Table 3, checked against Appendix C)

| Group_Model | Stress balance (Table 3) | Resolution km (Table 3) | Appendix C says | Discrepancy? |
|---|---|---|---|---|
| AWI_PISM1 | Hybrid | 8 | 8 km (16 km early in spin-up); SIA+SSA hybrid | – |
| DOE_MALI | HO | 2–20 | 3-D first-order Blatter–Pattyn; 2 km at grounding lines / marine areas, 20 km interior, ≤6 km on large shelves | – |
| ILTS_PIK_SICOPOLIS | Hybrid | 8 | Hybrid SIA/SStA; 8 km | – |
| IMAU_IMAUICE1 | Hybrid | 32 | SIA+SSA; 32 km | – |
| IMAU_IMAUICE2 | Hybrid | 32 | same | – |
| JPL1_ISSM | SSA | 2–50 | 2-D SSA; **1 km** at coast to 50 km interior; ≤8 km on shelves | **Yes:** 2 km (Table 3) vs 1 km (C5) |
| LSCE_GRISLI | Hybrid | 16 | SIA + SSA-as-sliding-law; 16 km | – |
| NCAR_CISM | L1L2 | 4 | Depth-integrated higher-order (Goldberg 2011); 4 km | – |
| PIK_PISM1 | Hybrid | 8 | SIA+SSA; 8 km | – |
| PIK_PISM2 | Hybrid | 8 | same | – |
| UCIJPL_ISSM | HO | 3–50 | Higher-order (Pattyn 2003); 3 km at margins to **30 km** inland | **Yes:** 50 (Table 3) vs 30 km (C9) |
| ULB_fETISh_16km | Hybrid | 16 | SSA for sliding + SIA for deformation; 16 km | – |
| ULB_fETISh_32km | Hybrid | 32 | C10 mentions only 16 km | 32 km variant not described in the text |
| UTAS_ElmerIce | Stokes | 4–40 | Full Stokes; ~4 km at grounding lines to ~40 km interior | – |
| VUB_AISMPALEO | SIA+SSA | 20 | SIA sheet + SSA shelf with a 1-cell transition zone (not a hybrid); 20 km | – |
| VUW_PISM | Hybrid | 16 | Spin-up refined 32 → 16 → 10 → **8 km** for the historical run | **Yes:** 16 (Table 3) vs 8 km (C13) |

---

## 4. Mapping and structure (both app and notebook)

1. **Substring fallback in `get_ism_meta`.** If no exact (Group, Model) key exists, it returns the first entry whose *group* string contains, or is contained in, the requested group, whatever the model and ice sheet. For AIS this fires only for ULB (#9), and lands on a plausible entry by luck. It's the mechanism that has mislabeled runs before, as the comment above the JPL/UCIJPL GIS keys says. **Recommend:** exact-key lookup only, and a startup `assert` that every run has an exact metadata row.
2. **No ice-sheet dimension in the key** (#10). **Recommend:** key on (ice_sheet, Group, Model), or, better, key on (publication, Group, Model).
3. **Two hand-maintained copies** (data.py and notebook cell 10) are identical today, but nothing keeps them in sync. **Recommend:** one metadata file (e.g. `metadata/ismip6_ais_models.csv` and `metadata/ismip6_ais_experiments.csv`) with a `source_citation` column (table / appendix section) per value. Both the app and the notebook would read it.
4. **Scenario is defined twice** (#13). `ismip6_ais_to_csv` writes a correct `RCP` column; `ais_exp_meta` holds a wrong `scenario`. **Recommend:** keep a single experiment table and drop the other.
5. **Missing experiment-level fields** that matter to the bias analysis: `melt_framework` (open/standard), `ocean_sensitivity` (Low/Medium/High/PIGL), `ice_shelf_fracture` (bool), `tier`. **Recommend** adding all four. Also make melt framework available as a grouping/ANOVA factor, or at least a filter.
6. **Melt parameterization depends on model × framework**, not on the experiment alone (Table 3's last two columns). E.g. PIK_PISM's open runs use PICO and UCIJPL's use PICOP; standard runs are "Nonlocal", "Local anom." etc. by model. **Recommend:** a field `melt_parameterization` computed from (model, framework).

---

## 5. Places where the paper disagrees with itself

These need a decision, or a cross-check against the model's own description paper.

| Item | Table 3 | Elsewhere in the paper |
|---|---|---|
| VUW_PISM initial year | 2015 | Table 4 lists a `historical` run; C13 describes a GCM-forced historical run 1950–2000. The archived series starts in 2001. |
| VUW_PISM resolution | 16 km | C13: final grid 8 km |
| IMAU_IMAUICE2 initial year | 1978 | C4: steady state assigned to **1900**, then a 79-yr forced run before 1979. (The archived series starts in 1980 for both IMAU runs.) |
| JPL1_ISSM init code | DA | C5 describes a 2-yr relaxation after the inversion, which is the caption's definition of **DA+** |
| UTAS_ElmerIce init code | DA | C11 describes a 10-yr relaxation after the inversion (again DA+ by definition) |
| Eq vs SP | AWI = Eq; PIK = SP; VUB = SP | C8 calls PIK's an "equilibrium simulation", and C12 calls VUB's a run "until a steady state is reached". Both read like Eq. |
| JPL1_ISSM resolution | 2–50 km | C5: 1 km at coast |
| UCIJPL_ISSM resolution | 3–50 km | C9: 3–30 km |
| UCIJPL_ISSM initial year | 2007 | Archived series starts in 2005 (minor; probably the archive's own hist start) |

---

## 6. Categorization decisions needed

These are the judgment calls. Each choice changes which runs are grouped together in the "what drives bias" ANOVA. All examples are from this paper.

### 6.1 Initialization: **DECIDED, Option C** (see the Decisions log and the class mapping below)

Table 3 uses six codes. The app collapses them to two (DA / Spin-up) and currently puts some runs on the wrong side even by its own scheme.

| Code (Table 3) | Meaning (caption) | AIS models |
|---|---|---|
| SP | spin-up | IMAU_IMAUICE2, PIK_PISM1, PIK_PISM2, VUB_AISMPALEO, VUW_PISM |
| Eq | equilibrium state | AWI_PISM1, IMAU_IMAUICE1 |
| SP+ | spin-up with ice-thickness targets (Pollard & DeConto 2012a) | ILTS_PIK_SICOPOLIS, LSCE_GRISLI, NCAR_CISM |
| DA\* | data assimilation of geometry only | ULB_fETISh_16km, ULB_fETISh_32km |
| DA+ | data assimilation + relaxation | DOE_MALI (and arguably JPL1, UTAS, per §5) |
| DA | data assimilation | JPL1_ISSM, UCIJPL_ISSM, UTAS_ElmerIce |

Why this is hard:
- **ULB's "DA\*" and NCAR/ILTS's "SP+" are essentially the same method.** Both nudge to observed thickness the Pollard & DeConto (2012a) way: C10 cites that procedure for ULB, and the Table 3 caption cites it for SP+. Yet the codes put them on opposite sides of a two-class split. LSCE's iterative drag inversion to match thickness (Le Clec'h et al. 2019) belongs to the same family. NCAR's own text calls it a hybrid.
- **Eq vs SP looks like an author choice, not a methodological difference** (see the §5 row).
- Seroussi et al. themselves contrast "steady state or spin-up" models with "data assimilation" models when discussing control-run drift (Sect. 4.2), which supports a split along those lines.

Options:
- **A. Keep the six paper codes verbatim.** Most faithful, but thin cells: Eq has 2 models and DA\* is effectively 1 model at 2 resolutions.
- **B. Two classes, following the paper's codes:** Spin-up = {SP, Eq, SP+}; Data assimilation = {DA, DA+, DA\*}. This is what Table 3 implies. It puts ULB on the DA side and NCAR/LSCE/ILTS on the spin-up side, even though their methods are nearly identical.
- **C. Three classes, by method:** (i) spin-up / equilibrium = {SP, Eq}; (ii) **thickness-nudged hybrid** = {SP+, DA\*} → ILTS, LSCE, NCAR, ULB×2; (iii) **velocity-based inversion** = {DA, DA+} → DOE, JPL1, UCIJPL, UTAS.
- Under any option, I'd store the verbatim paper code **as well as** the class, so the class can be re-derived later.

**Class mapping for ISMIP6-AIS under Option C**

| Class (`init_class`) | Paper codes | AIS models | Runs |
|---|---|---|---|
| Spin-up / equilibrium | SP, Eq | AWI_PISM1, IMAU_IMAUICE1, IMAU_IMAUICE2, PIK_PISM1, PIK_PISM2, VUB_AISMPALEO, VUW_PISM | 63 |
| Hybrid (thickness-nudged) | SP+, DA\* | ILTS_PIK_SICOPOLIS, LSCE_GRISLI, NCAR_CISM, ULB_fETISh_16km, ULB_fETISh_32km | 85 |
| Data assimilation (velocity inversion) | DA, DA+ | DOE_MALI, JPL1_ISSM, UCIJPL_ISSM, UTAS_ElmerIce | 40 |

Notes on the edge cases:
- **VUW_PISM stays in Spin-up.** Its multistage spin-up resets ice thickness to observed values at each grid refinement (C13). That pulls the geometry toward observations, but it doesn't tune basal friction or melt to fit them, and drift is removed by the resets. This matches the SP code.
- **ILTS_PIK_SICOPOLIS is a Hybrid**, and arguably the most data-constrained of them. Besides the topography nudging, it tunes basal sliding coefficients per IMBIE basin against observed velocities (C3).
- **ULB's runs (42) are now in a different class than before** (Data assimilation → Hybrid), and so are NCAR and LSCE. The Hybrid class holds 85 of 188 AIS runs, 42 of them ULB. This makes the model-weighting question in §6.5 more pressing.
- Compared with the current app labels, **7 of 16 AIS models change class**: IMAU_IMAUICE1/2 (Data assimilation → Spin-up); LSCE_GRISLI, NCAR_CISM and ULB ×2 (Data assimilation → Hybrid); ILTS_PIK_SICOPOLIS (Spin-up → Hybrid). The other 9 stay on the same side and only get the new class names. That relabels 105 of 188 runs (20 + 12 + 19 + 42 + 12).

The choice also has to work for Goelzer 2020 GIS, which uses its own codes, and for the external studies. The class scheme you pick here will be applied to every paper.

### 6.2 Sliding law: **DECIDED, four-class cascade** (see the Decisions log; literature check in §9)

What Seroussi 2020 actually states:

| Paper description | Models |
|---|---|
| Linear | DOE_MALI, UTAS_ElmerIce |
| Power law / Weertman | LSCE_GRISLI, ULB_fETISh ×2 ("coefficient of 2") |
| Weertman with effective-pressure dependence (Budd-type / height above buoyancy) | ILTS_PIK_SICOPOLIS, VUB_AISMPALEO |
| Pseudo-plastic: power law + Mohr–Coulomb till yield stress | PIK_PISM1/2 |
| Regularized Coulomb (Schoof 2005) | NCAR_CISM |
| **Not stated** | AWI_PISM1, IMAU_IMAUICE1/2, JPL1_ISSM, UCIJPL_ISSM, VUW_PISM (74 of 188 runs) |

Options:
- **A. By exponent / functional form:** Linear / Nonlinear power law (Weertman) / Coulomb-limited (pseudo-plastic + regularized Coulomb). Three classes.
- **B. By whether basal stress saturates** (the property most often linked to marine ice-sheet retreat sensitivity): *Power law incl. linear* vs *Coulomb-limited*. Two classes; DOE, UTAS, LSCE, ULB, ILTS and VUB on one side, PIK and NCAR on the other.
- **C. Finer:** A, plus a separate yes/no field for effective-pressure dependence (ILTS, VUB, PIK, NCAR yes).
- **"Not stated" runs:** either (i) label them "Not stated in source" and leave them out of the sliding-law ANOVA, or (ii) fill them from the cited model-description papers, and mark them as secondary-source values with a citation. Candidates: AWI_PISM → Kleiner & Humbert / PISM docs; JPL1 → Schlegel et al. 2018; UCIJPL → Morlighem et al. 2013; VUW → Golledge et al. 2019; IMAU → de Boer et al. 2014. Option (ii) means reading those papers.

### 6.3 Stress balance: **DECIDED, Option C (Shallow / Higher-order)**

The Table 3 codes give six categories: Hybrid (9 models), SIA+SSA (VUB), SSA (JPL1), L1L2 (NCAR), HO (DOE, UCIJPL), Stokes (UTAS). Four of them have one or two models.

Options:
- **A. Verbatim:** six classes, several singletons.
- **B. Grouped by order:** (i) shallow hybrids = Hybrid + SIA+SSA; (ii) SSA-only = JPL1; (iii) depth-integrated higher-order = L1L2; (iv) 3-D higher-order = HO + Stokes.
- **C. Two classes:** shallow (Hybrid, SIA+SSA, SSA) vs higher-order (L1L2, HO, Stokes).

One sub-question: should VUB's "SIA+SSA" count as a hybrid? Its grounded ice is SIA-only with Weertman sliding, which is physically different from PISM-style hybrids where the SSA itself is the sliding solution. The paper keeps the two separate.

**Class mapping for ISMIP6-AIS under Option C.** Each model was cross-checked against S19, ABU, LAR and Sch18 (sources in §9a).

| Class (`stress_balance_class`) | Model: paper code (other sources) | Runs |
|---|---|---|
| Shallow | AWI_PISM1: Hybrid (LAR: hybrid) · ILTS_PIK_SICOPOLIS: Hybrid (S19, LAR, ABU: hybrid) · IMAU_IMAUICE1/2: Hybrid (S19, LAR: SIA+SSA; ABU: Hybrid) · JPL1_ISSM: SSA (S19, LAR, ABU: SSA; Sch18: SSA for forward runs, L1L2 only during initialization) · LSCE_GRISLI: Hybrid (S19, LAR, ABU) · PIK_PISM1/2: Hybrid (S19, LAR, ABU) · ULB_fETISh_16/32km: Hybrid (S19, LAR, ABU) · VUB_AISMPALEO: SIA+SSA (S19, LAR: SIA sheet + SSA shelf) · VUW_PISM: Hybrid (LAR) | 141 (12 configurations) |
| Higher-order | DOE_MALI: HO (S19, LAR: Blatter–Pattyn; ABU: LMLa) · NCAR_CISM: L1L2 (S19, LAR: depth-integrated HO; ABU: L1L2) · UCIJPL_ISSM: HO (S19, LAR: Pattyn 2003) · UTAS_ElmerIce: Stokes (S20 only; UTAS isn't in S19/ABU/LAR) | 47 (4 configurations) |

All sources agree on every model's stress balance; no corrections were needed beyond adding the field.

Notes:
- **JPL1_ISSM (SSA) counts as Shallow.** SSA has no vertical shear in grounded ice, so it is shallower than the Hybrid models in that respect. Option C puts it with them.
- **NCAR_CISM (L1L2) counts as Higher-order.** L1L2 is depth-integrated but keeps the vertical-shear terms.
- **Collinearity with initialization:** 3 of the 4 Higher-order configurations (DOE, UCIJPL, UTAS) are also the velocity-inversion data-assimilation models. The fourth (NCAR) is a Hybrid. The only velocity-inversion model that is Shallow is JPL1. So in the AIS ensemble, stress balance and initialization class largely move together, and an ANOVA will have trouble attributing bias to one rather than the other.

### 6.4 Resolution: **DECIDED, Option B (binned `res_min_km`)**

Eleven models use uniform grids (4–32 km). Five use variable meshes (DOE 2–20, JPL1 2–50, UCIJPL 3–50, UTAS 4–40; in Table 3 the ranges are fine near grounding lines and coarse inland).

Options:
- **A.** Store `res_min_km` and `res_max_km`, and analyze with `res_min_km` (the finest resolution, at grounding lines, which is what controls grounding-line migration). Treat it as a continuous (log) covariate in the regression.
- **B.** The same values, but binned for grouping/ANOVA, e.g. ≤5, 5–10, 10–20, >20 km. With `res_min`: ≤5 = DOE, JPL1, NCAR, UCIJPL, UTAS; 5–10 = AWI, ILTS, PIK×2; 10–20 = LSCE, ULB16, VUB, VUW; >20 = IMAU×2, ULB32.
- For the §5 discrepancies (JPL1 1 vs 2 km, VUW 16 vs 8 km): use Table 3, per the agreed ranking, and note the alternative.

**Bin mapping for ISMIP6-AIS** (`res_min_km` from Table 3, per the agreed ranking; cross-checked against S19, ABU, LAR and Sch18, sources in §9a)

| `res_bin` | Model (res_min_km) | Runs |
|---|---|---|
| ≤5 km | DOE_MALI (2) · JPL1_ISSM (2) · NCAR_CISM (4) · UCIJPL_ISSM (3) · UTAS_ElmerIce (4) | 59 (5 configurations) |
| >5–10 km | AWI_PISM1 (8) · ILTS_PIK_SICOPOLIS (8) · PIK_PISM1 (8) · PIK_PISM2 (8) | 41 (4) |
| >10–20 km | LSCE_GRISLI (16) · ULB_fETISh_16km (16) · VUB_AISMPALEO (20) · VUW_PISM (16) | 47 (4) |
| >20 km | IMAU_IMAUICE1 (32) · IMAU_IMAUICE2 (32) · ULB_fETISh_32km (32) | 41 (3) |

**Cross-check results**
- **Agree everywhere:** DOE (2–20), ILTS (8), IMAU (32), LSCE (16), NCAR (4), VUB (20), ULB_16km (16), UTAS (4–40; S20 only).
- **JPL1_ISSM:** Table 3 gives 2 km. S20 C5, S19 and LAR all say **1 km** at the coast, and ABU gives 1–50. The bin is ≤5 km either way, but `res_min_km` should arguably be 1. I'm keeping 2 per the tables-first rule and noting the alternative.
- **UCIJPL_ISSM:** Table 3 gives 3–50. S20 C9 and S19 say 3–**30**, and LAR says 3–50. `res_min_km` = 3 either way; only `res_max_km` is in question.
- **VUW_PISM: this changes the bin.** Table 3 gives **16 km**, but S20 C13 says the grid was refined to **8 km** for the historical run that leads into the projections. The other sources use different VUW setups (LAR 16 km with PISM 0.7.1; G19 5 km). With 16 the model is in >10–20 km; with 8 it is in >5–10 km. Tables-first gives 16. This needs your call (§5) and affects 4 runs.
- **PIK_PISM:** S20 says 8 km. ABU and LAR describe a different 4 km PIK setup, so they don't apply here.
- **AWI_PISM1:** S20 says 8 km (16 km only early in the spin-up). LAR's 16 km is AWI's LARMIP setup, not ISMIP6.
- **ULB_fETISh_32km:** only Table 3 mentions the 32 km variant. The data archive's model name (`fETISh_32km`) independently confirms it.

**Collinearity (important).** The ≤5 km bin is almost exactly the same set of models as the Higher-order class (DOE, NCAR, UCIJPL, UTAS) plus JPL1. It also overlaps heavily with the velocity-inversion initialization class (DOE, JPL1, UCIJPL, UTAS). In the AIS ensemble, finest resolution, stress balance and initialization are three nearly interchangeable descriptions of the same five models. An ANOVA with all three as factors can't attribute bias among them. Consider fitting them one at a time, or reporting their combined effect.

### 6.5 Which ISMIP6-AIS runs belong in the bias analysis at all

This isn't strictly metadata, but it depends on the new fields:
- ~~**Ice-shelf collapse** runs (exp11/12; 14 runs)~~: **DECIDED, excluded from the bias analysis.** Runs excluded: AWI_PISM1, UCIJPL_ISSM, ULB ×2 (exp11 and exp12); DOE_MALI, ILTS_PIK_SICOPOLIS, IMAU ×2, JPL1_ISSM, LSCE_GRISLI (exp12 only). NCAR_CISM, PIK ×2, UTAS, VUB and VUW ran neither.
- ~~**Ocean-sensitivity variants** (exp09/10/13; 37 runs)~~: **DECIDED, kept.** Side effect: all three are NorESM1-M RCP8.5 standard-framework runs, the same forcing as exp05. So NorESM1-M RCP8.5 accounts for 4 of the standard experiments per model, and it gets more weight than other GCMs in the `climate_model` factor.
- ~~**Melt framework**~~: **DECIDED, an ANOVA factor (`ocean_forcing`), with `ocean_sensitivity` as a second factor.** Levels in the AIS bias set (174 runs):

  | | High | Low | Medium | PIGL | Total |
  |---|---|---|---|---|---|
  | Open | 0 | 0 | 48 | 0 | 48 |
  | Standard | 12 | 12 | 89 | 13 | 126 |

  How to read these factors:
  - **`ocean_sensitivity` is nested in Standard forcing.** High, Low and PIGL occur only in Standard runs, and only with NorESM1-M RCP8.5 (exp09/10/13). Their contrast with Medium is effectively within-model, within-GCM, which is clean. But their effect can't be estimated for Open runs or for other GCMs.
  - **`ocean_forcing` is partly confounded with model.** Only 5 configurations ran both frameworks (AWI_PISM1, NCAR_CISM, UCIJPL_ISSM, ULB ×2). Three ran Open only (PIK_PISM1, PIK_PISM2, VUW_PISM), all of them PISM, Coulomb-limited, Spin-up. Eight ran Standard only. So outside those 5 models, an Open-vs-Standard difference is partly a difference between PISM and everything else.
  - **"Open" isn't one melt scheme.** Each group used its own (Table 3: Quad, PICO, PICOP, Plume, Lin, Nonlocal + Slope). The per-model scheme is stored as `melt_parameterization` (§4.6), so it can be inspected, but it is not a factor.
- ~~**Model weighting**~~: **DECIDED, per run (unchanged).** ULB ×2 = 38 of 174 bias-set runs (22 %), AWI_PISM1 = 19 (11 %), and UTAS_ElmerIce = 3 (2 %). Results for the classes these models sit in (ULB is all of "Power law") are weighted toward them.

---

## 7. Recommended changes

### 7a. Unambiguous now (no category decision needed)

1. `ais_exp_meta` scenario: expA1, A2, A3, A5, A6, A7 → **RCP8.5**; expA4, A8 → **RCP2.6**.
2. `ais_exp_meta` climate_model: "NorESM" → **"NorESM1-M"**.
3. Replace `basal_melt_param` with `ocean_sensitivity`: exp09 **High**, exp10 **Low**, exp13 **PIGL**, all other runs (open included, per Table 1) **Medium**. Add `ocean_forcing` (Open/Standard; replaces `protocol`), `ice_shelf_fracture` (True for exp11, exp12) and `tier`. Add `ocean_forcing` and `ocean_sensitivity` to `GROUP_DIMENSIONS` and `ANOVA_CHARACTERISTICS`.
4. Sliding law where Appendix C is explicit: DOE_MALI → **Linear**; UTAS_ElmerIce → **Linear**; ULB_fETISh → **Power law (Weertman)**. Final class labels follow decision 6.2.
5. Add exact keys `("ULB","fETISh_16km")` and `("ULB","fETISh_32km")`.
6. Remove the substring fallback; assert an exact match for every run.
7. Key model metadata by ice sheet (or publication), so the AIS fixes don't change GIS runs.
7b. In `analysis.run_weights`, compute `typical` from the number of distinct (Group, Model) configurations, not distinct groups (`len(set(_GROUP[ism]))` → count of unique `group`+`model` pairs). This follows the configuration-independence decision.
8. Add `stress_balance` (verbatim Table 3 code), `res_min_km`, `res_max_km`, `init_code` (verbatim Table 3) and `init_year` (Table 3) from §3b.
9. Fix the citation comment: "Seroussi et al. 2020 Table A1" → "Table 3 and Appendix C".
10. Make the same edits in the notebook, or better, have the notebook and app read one shared metadata file (§4.3).

### 7b. Pending your decisions

- ~~`initialization` class (6.1)~~: **decided (Option C)**. Classes are now in §7c.
- `sliding_law` class and how to handle the 6 "not stated" models (6.2).
- ~~`stress_balance` grouping (6.3)~~: **decided (Option C)**. ~~Resolution (6.4)~~: **decided (Option B, binned `res_min_km`)**.
- ~~Run inclusion and weighting (6.5)~~: **all decided** (collapse runs excluded; sensitivity variants kept; ocean forcing and sensitivity as factors; per-run weights).
- Resolving the internal paper inconsistencies in §5.

### 7c. Proposed AIS model record (verbatim-paper values; `init_class` and `stress_balance_class` filled; `res_bin` filled; sliding-law fields in §9d)

| Group | Model | ice_model | stress_balance | stress_balance_class | res_min_km | res_max_km | res_bin | init_code | init_year | init_class | sliding_law_paper | open_melt | standard_melt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AWI | PISM1 | PISM | Hybrid | Shallow | 8 | 8 | >5–10 | Eq | 2005 | Spin-up | Not stated | Quad | Nonlocal |
| DOE | MALI | MALI | HO | Higher-order | 2 | 20 | ≤5 | DA+ | 2015 | DA (velocity) | Linear | N/A | Nonlocal anom. |
| ILTS_PIK | SICOPOLIS | SICOPOLIS | Hybrid | Shallow | 8 | 8 | >5–10 | SP+ | 1990 | Hybrid | Weertman–Budd, sub-melt sliding | N/A | Nonlocal |
| IMAU | IMAUICE1 | IMAU-ICE | Hybrid | Shallow | 32 | 32 | >20 | Eq | 1978 | Spin-up | Not stated | N/A | Local anom. |
| IMAU | IMAUICE2 | IMAU-ICE | Hybrid | Shallow | 32 | 32 | >20 | SP | 1978 | Spin-up | Not stated | N/A | Local anom. |
| JPL1 | ISSM | ISSM | SSA | Shallow | 2 | 50 | ≤5 | DA | 2007 | DA (velocity) | Not stated | N/A | Nonlocal |
| LSCE | GRISLI | GRISLI | Hybrid | Shallow | 16 | 16 | >10–20 | SP+ | 1995 | Hybrid | Power law (Weertman) | N/A | Nonlocal |
| NCAR | CISM | CISM | L1L2 | Higher-order | 4 | 4 | ≤5 | SP+ | 1995 | Hybrid | Regularized Coulomb (Schoof 2005) | Nonlocal + Slope | Nonlocal |
| PIK | PISM1 | PISM | Hybrid | Shallow | 8 | 8 | >5–10 | SP | 1850 | Spin-up | Pseudo-plastic (Mohr–Coulomb till) | PICO | N/A |
| PIK | PISM2 | PISM | Hybrid | Shallow | 8 | 8 | >5–10 | SP | 2015 | Spin-up | Pseudo-plastic (Mohr–Coulomb till) | PICO | N/A |
| UCIJPL | ISSM | ISSM | HO | Higher-order | 3 | 50 | ≤5 | DA | 2007 | DA (velocity) | Not stated | PICOP | Nonlocal |
| ULB | fETISh_16km | f.ETISh (Kori-ULB) | Hybrid | Shallow | 16 | 16 | >10–20 | DA\* | 2005 | Hybrid | Power law (coefficient 2) | Plume | Nonlocal |
| ULB | fETISh_32km | f.ETISh (Kori-ULB) | Hybrid | Shallow | 32 | 32 | >20 | DA\* | 2005 | Hybrid | Power law (coefficient 2) | Plume | Nonlocal |
| UTAS | ElmerIce | Elmer/Ice | Stokes | Higher-order | 4 | 40 | ≤5 | DA | 2015 | DA (velocity) | Linear | N/A | Local |
| VUB | AISMPALEO | AISMPALEO | SIA+SSA | Shallow | 20 | 20 | >10–20 | SP | 2000 | Spin-up | Weertman (∝ 1/height above buoyancy) | N/A | Nonlocal anom. |
| VUW | PISM | PISM | Hybrid | Shallow | 16 | 16 | >10–20 | SP | 2015 | Spin-up | Not stated | Lin | N/A |

Source for every cell: Seroussi et al. (2020) Table 3, except `sliding_law_paper` (Appendix C1–C13). §5 lists the cells where the text disagrees.

---

## 8. Open items for the next pass

- **Configuration independence vs. earlier collapsing choices (raise at those papers):** `data.py` notes that Goelzer 2025's NORCE CISM resolution/tuning variants were collapsed to a single Model = "CISM" ("per the user's explicit choice"), and that Rahlves 2025's ERA5- and ESM-initialized branches share one metadata row. Both conflict with the decision that every configuration is independent, and collapsing the NORCE variants would also hide their different resolutions. Revisit both when reviewing those papers.

- Decide 6.1–6.5. The answers set the class columns and the code changes.
- If you choose option (ii) in 6.2, read the model-description papers for the six "not stated" sliding laws.
- Next paper in sequence: Goelzer et al. (2020), ISMIP6-Greenland. That review will also settle the shared AIS/GIS keys (#10) and the GIS substring-fallback cases I noticed in passing: `ILTS_PIK/SICOPOLIS1`, `SICOPOLIS2` and `AWI/ISSM3` have no exact keys.

---

## 9. Sliding laws: literature check (2026-11-23)

As requested, every AIS model's sliding law was checked beyond Seroussi 2020: the six "not stated" models were filled from other papers, and the ten stated ones were cross-checked against other papers using the same model from the same group. Classes follow the agreed cascade: (1) Coulomb-limited → (2) Budd (depends on effective pressure or height above flotation) → (3) Linear → (4) Power law.

### 9a. Sources used

| Key | Paper | What it gives |
|---|---|---|
| S20 | Seroussi et al. 2020, TC 14, 3033 (Appendix C) | The ISMIP6 configurations themselves |
| S19 | Seroussi et al. 2019, initMIP-Antarctica, TC 13, 1441 (Appendix B) | The same groups one experiment earlier. Sliding-law text is nearly identical to S20 for most groups. |
| ABU | Sun et al. 2020, ABUMIP, J. Glaciol. 66, 891 (Table 2, Appendix A) | **The only source that tabulates m and q per model**, and contemporaneous with ISMIP6 |
| LAR | Levermann et al. 2020, LARMIP-2, ESD 11, 35 (Appendix A1–A16) | Per-model descriptions from the same groups, often in different words than S20 |
| dB15 | de Boer et al. 2015, PLISMIP-ANT, TC 9, 881 (Table 3) | The IMAU (ANICE) sliding law |
| Q18 | Quiquet et al. 2018, GRISLI 2.0, GMD 11, 5003 (Sect. 2.1.3) | The GRISLI friction equation |
| LC19 | Le clec'h et al. 2019, GMD 12, 2481 (Sect. 2) | The GRISLI inversion and its linear viscous till |
| Sch18 | Schlegel et al. 2018, TC 12, 3511 (Sect. 2.1) | The JPL ISSM friction equation (S20 cites it for JPL1_ISSM) |
| G19 | Golledge et al. 2019, Nature 566, 65 (Methods, accepted MS) | The VUW PISM setup |
| GW | Gladstone & Wang, Elmer/Ice Antarctic inversion methods note (Zenodo 5862046) | Elmer/Ice sliding relation (Gladstone is a UTAS co-contributor) |

Not available: AWI's EGU2020 slides on its ISMIP6 PISM setup (file over the fetch size limit); the Golledge 2015 Nature Methods (paywalled); Morlighem et al. 2013 JGR (paywalled).

### 9b. Resolved sliding law per model

| Model | S20 says | Other sources | Resolved law | **Class** | Confidence | App now | Change? |
|---|---|---|---|---|---|---|---|
| AWI_PISM1 | Not stated | LAR A10 (AWI, equilibrium spin-up): pseudo-plastic friction law, depending on till yield strength and stored basal water (Bueler & Brown 2009; Winkelmann et al. 2011). ABU: AWI's *paleo* variant, Coulomb q = 0.6. S19 B2: not stated. | Pseudo-plastic, Mohr–Coulomb till; q not known for this configuration | **Coulomb-limited** | Medium. Same group and model, but no source for the exact ISMIP6 8 km PISM 1.1.4 setup. | Pseudo-plastic | Label OK |
| DOE_MALI | Linear | S19 B5: linear. LAR A9: linear. ABU: Weertman m = 1. | Linear, spatially varying coefficient, no N | **Linear** | High (4 sources agree) | **Coulomb** | **Fix** |
| ILTS_PIK_SICOPOLIS | Weertman–Budd, sub-melt sliding | S19 B7: Weertman-type with sub-melt sliding. LAR A15: Weertman–Budd. ABU: m = 3, p = 2 (the only model ABU names with N dependence). | Weertman–Budd, m = 3, p = 2 | **Budd** | High | Weertman | **Class changes** |
| IMAU_IMAUICE1/2 | Not stated | S19 B8 and LAR A6: not stated (SSA velocities added as sliding). ABU Table 2: **Coulomb q = 0**, and the text calls IMAU-ICE the only model with a Coulomb plastic sliding law. dB15 Table 3, ANICE (the SIA+SSA model): **"Mohr–Coulomb plastic law"**. | Plastic Mohr–Coulomb (q = 0) | **Coulomb-limited** | Medium-high (2 independent same-group sources agree) | **Weertman** | **Fix**. The "Weertman" label probably comes from dB15's *SIA-only* variant, which that paper also calls "IMAU-ICE". |
| JPL1_ISSM | Not stated | **Sch18 (cited by S20 for this setup): Budd-type law, τ_b = −α² N v_b, with N = g(ρ_i H + ρ_w z_b).** LAR A7: Budd (Budd et al. 1979). ABU: Weertman m = 1 (ABU says most models ignore N; Sch18 contradicts this for JPL). | Linear Budd, N from ice overburden and bed depth | **Budd** | High (primary model paper) | Budd / Schoof | **Drop "Schoof"** |
| LSCE_GRISLI | Power law (Weertman 1957) | S19 B10 and LAR A5: same wording. ABU: **Coulomb q = 1** (i.e. linear). **Q18: τ_b^m = −β u_b with m = 1 ("viscous deformation") and β = C_f·N.** **LC19: "linear viscous till", τ_b = −β u_b, β varies in space and is inverted.** S20 C6: β inverted, then held constant. | Nominal power law with m = 1 (so linear). In ISMIP6, β was inverted and frozen, so N does not act during the run. | **Linear** | Medium. If the ISMIP6 runs kept β = C_f·N with N evolving, the class would be **Budd**. Neither S20 nor LC19 says so. | **Weertman** | **Fix** |
| NCAR_CISM | Like Schoof 2005 | S19 B11 and LAR A3: same. ABU: Weertman m = 3 / Coulomb. | Regularized Coulomb (Schoof 2005), m = 3 | **Coulomb-limited** | High | Regularized Coulomb | OK |
| PIK_PISM1/2 | Power law + Mohr–Coulomb till yield stress | S19 B12: Mohr–Coulomb. ABU: pseudo-plastic q = 0.75. LAR A12 (a different 4 km setup): "plastic", independent of velocity. | Pseudo-plastic Mohr–Coulomb (q = 0.75 in ABU) | **Coulomb-limited** | High for the class; q depends on the setup | Pseudo-plastic | OK |
| UCIJPL_ISSM | Not stated | S19 B14: not stated. **LAR A8 (UCI): inversion "based on a Budd friction law (Budd et al., 1979)".** Same ISSM code and the same N form as Sch18. | Linear Budd | **Budd** | Medium-high (one same-group source; Morlighem 2013 not accessible) | **Regularized Coulomb** | **Fix** |
| ULB_fETISh_16/32km | Power law, coefficient 2 | S19 B15 and LAR A4: same. ABU: Weertman m = 2. ABU also notes ULB re-ran ABUK with m = 1–4 as a sensitivity test, not in ISMIP6. | Weertman, m = 2 | **Power law** | High | Weertman / Coulomb | **Drop "Coulomb"** |
| UTAS_ElmerIce | Linear | Not in S19, ABU or LAR. GW: whole-Antarctic Elmer/Ice inversions use a linear sliding relation, starting from Gladstone et al. 2019 / S20. | Linear | **Linear** | Medium-high | **Regularized Coulomb** | **Fix** |
| VUB_AISMPALEO | Weertman ∝ 1/height above buoyancy | S19 B16 and LAR A1: the same sentence. | Weertman with height-above-buoyancy dependence; exponent not stated | **Budd** | High for the form | Weertman | **Class changes** |
| VUW_PISM | Not stated (refers to Golledge 2019) | G19 Methods: PISM 0.7.3; enhancement factors only, no sliding law. LAR A13: not stated. S19 B1 (ARC = VUW's Antarctic Research Centre, Golledge): not stated. **ABU, ARC-PISM (Golledge): Coulomb q = 0.75.** | Pseudo-plastic q = 0.75 | **Coulomb-limited** | Medium (one same-PI source; nothing specific to the ISMIP6 runs) | Pseudo-plastic (Coulomb) | Label OK |

### 9c. Result under the cascade

| Class | Models | Runs (of 188) |
|---|---|---|
| Coulomb-limited | AWI_PISM1, IMAU_IMAUICE1, IMAU_IMAUICE2, NCAR_CISM, PIK_PISM1, PIK_PISM2, VUW_PISM | 72 |
| Budd | ILTS_PIK_SICOPOLIS, JPL1_ISSM, UCIJPL_ISSM, VUB_AISMPALEO | 51 |
| Power law | ULB_fETISh_16km, ULB_fETISh_32km | 42 |
| Linear | DOE_MALI, LSCE_GRISLI, UTAS_ElmerIce | 23 |

**Changes from the current app labels:** 7 of 16 configurations had a sliding law that was wrong or unsupported: DOE, IMAU ×2, LSCE, UCIJPL and UTAS were wrong, and JPL1's "/ Schoof" was unsupported. ULB's "/ Coulomb" was also unsupported. Under the new classes, ILTS and VUB move from Weertman to Budd.

**Things to know before using this as an ANOVA factor:**
1. **"Power law" contains only ULB** (one model at two resolutions, 42 runs). In the AIS ensemble, the sliding-law effect for this class can't be separated from the ULB model itself.
2. **Several classes are confounded with model family.** All PISM configurations are Coulomb-limited, and both ISSM configurations are Budd. Sliding law and ice-model code will be strongly collinear in a multi-factor ANOVA.
3. **"Budd" groups linear and nonlinear laws.** JPL1/UCIJPL are linear in velocity, while ILTS is m = 3. If exponent matters more than N dependence for your question, the stored `sliding_exponent` field lets you regroup.

### 9d. Proposed sliding-law fields (replacing `sliding_law`)

| Group_Model | sliding_law_class | sliding_law_detail | sliding_exponent | effective_pressure | sliding_law_source | confidence |
|---|---|---|---|---|---|---|
| AWI_PISM1 | Coulomb-limited | Pseudo-plastic, Mohr–Coulomb till | q unknown | yes (till) | LAR A10 | medium |
| DOE_MALI | Linear | Linear, inverted coefficient | m = 1 | no | S20 C2; S19; LAR; ABU | high |
| ILTS_PIK_SICOPOLIS | Budd | Weertman–Budd, sub-melt sliding | m = 3, p = 2 | yes | S20 C3; ABU | high |
| IMAU_IMAUICE1 | Coulomb-limited | Plastic Mohr–Coulomb | q = 0 | yes (till) | ABU; dB15 | medium-high |
| IMAU_IMAUICE2 | Coulomb-limited | Plastic Mohr–Coulomb | q = 0 | yes (till) | ABU; dB15 | medium-high |
| JPL1_ISSM | Budd | Budd, τ_b = −α² N v_b | linear | yes (N = g(ρ_i H + ρ_w z_b)) | Sch18; LAR A7 | high |
| LSCE_GRISLI | Linear | Viscous till, inverted β held fixed | m = 1 | no in ISMIP6 runs (β = C_f·N in GRISLI 2.0) | Q18; LC19; ABU | medium |
| NCAR_CISM | Coulomb-limited | Regularized Coulomb (Schoof 2005) | m = 3 | yes | S20 C7; ABU | high |
| PIK_PISM1 | Coulomb-limited | Pseudo-plastic, Mohr–Coulomb till | q = 0.75 | yes (till) | S20 C8; ABU | high |
| PIK_PISM2 | Coulomb-limited | Pseudo-plastic, Mohr–Coulomb till | q = 0.75 | yes (till) | S20 C8; ABU | high |
| UCIJPL_ISSM | Budd | Budd (inverted) | linear | yes | LAR A8 | medium-high |
| ULB_fETISh_16km | Power law | Weertman | m = 2 | no | S20 C10; ABU | high |
| ULB_fETISh_32km | Power law | Weertman | m = 2 | no | S20 C10; ABU | high |
| UTAS_ElmerIce | Linear | Linear, inverted coefficient | m = 1 | no | S20 C11; GW | medium-high |
| VUB_AISMPALEO | Budd | Weertman ∝ 1/height above buoyancy | not stated | yes (height above buoyancy) | S20 C12; S19; LAR | high |
| VUW_PISM | Coulomb-limited | Pseudo-plastic, Mohr–Coulomb till | q = 0.75 | yes (till) | ABU (ARC-PISM) | medium |

---

## 10. Implementation (2026-11-23)

Applied in this repository on git branch **`metadata-review-ais`** (branched from `main` at `61898ab`), commit `7c99df3`; the per-study files of §10b followed in `2cf222f`.

| File | Change |
|---|---|
| `dash_app/metadata/ismip6_ais_models.csv` (new) | The 16 reviewed AIS configurations from §7c/§9d: verbatim Table 3 codes plus the decided classes (`initialization`, `sliding_law`, `stress_balance`, `resolution`) and sliding-law detail, exponent, effective pressure, source and confidence |
| `dash_app/metadata/ismip6_ais_experiments.csv` (new) | Table 1: CMIP5 RCP scenarios, NorESM1-M, `ocean_forcing`, `ocean_sensitivity` (open = Medium), `ice_shelf_fracture`, `tier` |
| `dash_app/metadata/unreviewed_models.csv`, `ismip6_gis_experiments_unreviewed.csv` (new) | ISMIP6-GIS and the external studies, carried over unchanged and keyed explicitly by (ice_sheet, publication, group, model). The 3 GIS rows that used to resolve through the substring fallback (AWI/ISSM3, ILTS_PIK/SICOPOLIS1, ILTS_PIK/SICOPOLIS2) keep their old values, with a note saying so. |
| `dash_app/metadata.py` (new) | The single loader shared by the app and the notebook. Exact (ice_sheet, group, model) lookups, with KeyError instead of a fallback. Unreviewed labels get an " (unreviewed)" suffix and new fields read "Not yet reviewed"; emulators and assessments read "Not applicable". |
| `dash_app/data.py` | Old lookup tables and `get_ism_meta`'s fallback removed. Experiment metadata keyed by (publication, exp). New run columns: `stress_balance`, `resolution`, `ocean_forcing`, `ocean_sensitivity`, `ice_shelf_fracture`, `review_status`. New `GROUP_DIMENSIONS`/`ANOVA_CHARACTERISTICS`, `CATEGORY_ORDER`, `BIAS_ELIGIBLE`. Hover shows classes, verbatim codes, resolution range, collapse flag and review status. |
| `dash_app/analysis.py` | `run_weights`: the typical ISMIP6 size is now runs per (group, model) configuration. `category_masks`: reviewed vocabularies in their natural order. |
| `dash_app/app.py` | Bias/ANOVA/low-bias section drops `~BIAS_ELIGIBLE` runs and says how many. Intro and weighting text updated. |
| `dash_app/tests/test_metadata.py` (new) | 14 tests pinning these decisions (all 41 tests pass) |
| `exploring_fits/analyze_slr_predictions_interactive.ipynb` | Cell 10 reads the shared metadata (the duplicate tables are gone). `get_ism_meta(group, model, ice_sheet)` is exact. Hover, `_sim_rows`, `GROUP_DIMENSIONS`, `INIT_COLOR` and `CHARACTERISTICS` carry the new fields. `misfit_dataframe` drops collapse runs. |
| `README.md` | Group-by list, metadata location, collapse exclusion, weighting wording |

**Checked in the running app:** AIS bias section shows "Based on 174 simulations … 14 runs with imposed ice-shelf collapse are left out", with all nine factors in the ANOVA. Initialization and sliding law come out as "confounded with others" in the joint model, as §6.3/§6.4 predicted. GIS runs unchanged on its unreviewed labels. No console errors.

**Checked in the notebook:** cells 1–31 ran without errors on a scratch copy (so no figure files were written into the repo). The AIS bias ANOVA uses 174 runs and all nine characteristics. One caution: in the notebook's *combined* statsmodels table, the AIS design is rank-deficient (residual df 151 = 174 − rank 23, though the terms list 28 df). That's the collinearity among model, sliding law, initialization, stress balance and resolution, and statsmodels still prints F values for terms the app labels "confounded with others". Read the one-way table, or the app's joint column, rather than those F values.

### 10b. Per-study metadata files (2026-11-23)

The metadata now lives in one `<Study>_models.csv` and one `<Study>_experiments.csv` per study in `dash_app/metadata/`, plus `encoding_rules.csv` (allowed values, display order and assignment rules) and `experiment_columns.csv` (how each experiments column is used, including verbatim-to-encoded value maps such as `Medium*` → `Medium`). In the models files, each characteristic has `_described` (source wording, quoted where possible), `_source`, `_encoded` and `_note` columns. `Seroussi2020_experiments.csv` is Table 1 verbatim. The nine other studies log their current encodings as they are (review_status "not yet reviewed" / "not applicable").

`metadata.py` validates every file on import. The run table is identical before and after the move in every analysis column, and the notebook's AIS ANOVA output is unchanged. Tests: 54 pass, including checks that every bundled run has a models and an experiments row, and that the bundled data's climate_model/scenario agree with each experiments file. Both READMEs now ask visitors to check the encodings and give the path.
