# Metadata review 11: Golledge et al. (2019), PISM Antarctica and Greenland (new study)

**Paper:** Golledge, N. R., Keller, E. D., Gomez, N., Naughten, K. A., Bernales, J., Trusel, L. D. and Edwards, T. L. (2019). Global environmental consequences of twenty-first-century ice-sheet melt. *Nature* 566, 65–72. https://doi.org/10.1038/s41586-019-0889-9. The accepted manuscript was read in full (https://nora.nerc.ac.uk/521866/1/SLR_FINAL.pdf). Citations to it below are "G19".

**Data:** the author's Open Science Framework project, https://osf.io/3edxf/ ("Ice sheet model outputs as described in Golledge et al., 2019, Nature 566, 65-72"; contributor N. Golledge). The paper doesn't link to it. Files used:
- `G19_Fig2-data.zip`: the Figure 2 series.
- `AIS-CTRL.nc` and `GIS-CTRL.nc`: the control runs.

**Licence:** **none stated** on OSF, which by default means all rights reserved. The user chose to add the data now (2026-11-24). Permission or a licence from the author should still be sought.

**Decision (user, 2026-11-24):** add now, restoring the control.

---

## 1. What was done and verified

| # | Item |
|---|---|
| 1 | **Runs.** Eight series, `sl_cont_<AIS\|GIS>-<45\|85>-<2100\|MWA-cdo>-bc.txt`: sea-level contribution in m, every 5 years 2000–2100. "MWA" is the meltwater-feedback runs: RCP8.5 MWA reaches 0.140 m (AIS) + 0.112 m (GIS) = 0.25 m in 2100, the paper's "up to 25 cm". "2100" is the standard-forcing runs. |
| 2 | **"bc" = control subtracted.** G19 Methods: "all of our sea-level-equivalent volume calculations are bias-corrected using a control run that also uses this initial state". |
| 3 | **Control restored.** See the note below the table. |
| 4 | **The control drifts a lot.** See the note below the table. |
| 5 | **Two published series start just off zero in 2000:** AIS RCP4.5 at −0.39 mm and GIS RCP4.5 at −0.15 mm. They're kept as published; the app zeroes every series at the window start. |
| 6 | **Interpolation.** The 5-yearly values are interpolated to annual, as for DeConto & Pollard 2016's decadal values. A 2015–2023 rate therefore rests on the 2015, 2020 and 2025 points. |
| 7 | **App.** "Golledge 2019 (AIS)" and "Golledge 2019 (GIS)", in publication order between DeConto & Pollard 2016 and Aschwanden 2019, coloured #aa88ee (CIEDE2000 ≥ 26 against every other study), unchecked by default. Each is weighted as one typical ISMIP6 configuration, since it isn't an ISMIP6 study. The notebook adds them as extra sources. |

**How the control is restored (item 3).**
- The control files hold 31 snapshots of `thk`, `topg` and `mask`, every 5 years from 1950 to 2100, matching the series' 5-year points.
- The control's ice mass above flotation comes from grounded cells (mask 2) with H − max(0, −topg) × 1028/910 > 0, times dx·dy × 910 kg/m³. These are PISM 0.7's default densities; the PISM command lines in the files set no others.
- Totals in 2000: 56.8 m sea-level equivalent (AIS) and 7.0 m (GIS).
- Restored mass = −(published m SLE) × 362.5 Gt/mm + (control − control in 2000). A test confirms the restored-minus-published difference is identical across each ice sheet's four runs.

**The control's drift (item 4).**
- The Antarctic control gains about 9.6 mm sea-level equivalent per decade (about 350 Gt/yr). The Greenland control loses about 1.4 mm per decade.
- With the control restored, the Antarctic runs gain 330–480 Gt/yr over 2015–2023, where IMBIE shows loss. The Greenland runs lose 330–650 Gt/yr.
- G19 notes for runs without the sub-grid melt scheme that "the control run also increases in volume". The control-subtracted series the paper reports remove this drift.

---

## 2. Model configuration (PISM 0.7.3, one row per ice sheet)

| Characteristic | Sources | Encoding |
|---|---|---|
| Ice model | "the Parallel Ice Sheet Model (PISM) version 0.7.3" | PISM |
| Stress balance | "a 'hybrid' ice sheet / shelf model that combines shallow approximations of the flow equations"; the control's command line reads `-stress_balance ssa` | **Shallow** |
| Sliding law | Control command lines: `-pseudo_plastic -pseudo_plastic_q 0.2` (AIS) or `0.6` (GIS), `-till_effective_fraction_overburden 0.02`, `-topg_to_phi`; Extended Data Fig. 5: "till friction angle in the piecewise linear basal traction parameter" | **Coulomb-limited** (cascade step 1, PISM Mohr–Coulomb till) |
| Initialization | See the note below the table. | **Hybrid (thickness-nudged)** (rule D1, as Seroussi 2020's VUW_PISM, whose spin-up follows this one); forward runs start 1860 |
| Resolution | "grid resolutions of 2.5 km (Greenland) and 5 km (Antarctica)"; control grid spacing 5.0 and 2.5 km | **≤5 km** (5 km is a bin edge and belongs to the lower bin) |

**Initialization.** A glacial-cycle spin-up with "hierarchical grid refinement" and "a novel 'nudging' technique, such that when we incrementally increase the resolution of the model domains we also reset bed topography and ice thicknesses to present-day values".

**Forcing:**
- **Climate model:** "the CMIP5 multi-model mean" for RCP4.5 and RCP8.5, encoded as "CMIP5 multi-model mean". The meltwater-feedback runs add LOVECLIM anomalies, shown in the hover.
- **Regional climate model:** "None (AOGCM output used directly)"; CMIP5 fields are applied on a RACMO2.3 climatology.
- **Ocean forcing:** **Open**, the group's own melt-factor parameterization.
- **Ocean sensitivity:** **Medium**, a single group calibration, as for DeConto 2021 (DC2).
- **Ice-shelf collapse:** none imposed.
