# Ice Sheet Simulation Explorer

**How do simulations of ice sheet mass loss compare with observations?**

This repository holds an interactive web app and an accompanying Jupyter notebook. Both compare published
simulations of the Antarctic and Greenland ice sheets with satellite observations of their mass loss. The
simulations come from ISMIP6 and several other ensembles; the observations are the IMBIE assessments.

The tool is meant for **user-motivated investigation**. Choose the studies, time windows and modeling
characteristics you're curious about, and see how each choice changes the comparison. The goal is to help
build intuition for, and understanding of, ice sheet modeling: which simulations reproduce recent observed
mass loss, what the best-matching simulations have in common, and how today's agreement (or disagreement)
relates to projections for 2100.

**Live app:** https://ice-sheet-simulation-explorer.onrender.com/
(on Render's free tier, so the first visit after a quiet spell can take a minute to wake up).

## The web app

The app (`dash_app/`) is a single scrolling page with shared controls in a sidebar:

- **Averaging window:** the years over which rates and bias are computed. Time series are zeroed at the
  window's start.
- **Simulation studies:** which ensembles to include (listed under [Data](#data)): the ISMIP6 ensembles
  first, then the other ice sheet model studies by publication date, then the Edwards et al. (2021) emulator
  results, which come from a statistical emulator rather than an ice sheet model.
- **Group by:** split every plot by study, ice sheet model, initialization, sliding law, climate
  scenario or GCM. You can optionally pool similar RCP and SSP scenarios into "composite" scenarios.
- **Units:** Gt of ice, or mm of sea-level equivalent.
- **Observations:** IMBIE3 (through 2023) or IMBIE2 (through 2020).

It has four sections:

1. **Rates of mass change:** the distribution of each simulation's average rate over the window, against
   the observed rate and its ±2σ uncertainty.
2. **Mass change through time:** median and 25–75% / 5–95% ranges of cumulative mass change, against
   the IMBIE observations. The plot covers the averaging window plus six years either side.
3. **Mass change by 2100:** the spread of projected 2015–2100 change for each group.
4. **What drives bias?** Bias is a simulation's rate minus the observed rate. This section uses a weighted
   ANOVA to show how much of the spread in bias each modeling characteristic explains. It also shows the
   bias by category, and which characteristics are over-represented among the best-matching 10% of
   simulations.

Large ensembles are weighted so that each study counts about as much as one ISMIP6 modeling group. That
keeps a 3,000-member emulator ensemble from swamping the dozen or so modeling groups behind each ISMIP6
ensemble.

### Comparing your own ensemble

Click **"+ Add your own ensemble"** in the sidebar to upload a CSV of your own simulations and see them
beside the published ones. The first column is the year, and each other column is one simulation holding
cumulative mass (or mass change) in Gt, with mass loss negative:

```
year,exp01,exp02,exp03
2015,0.0,0.0,0.0
2016,-231.4,-198.7,-260.2
...
```

After uploading, choose whether the simulations are of Antarctica or Greenland. Uploads stay in your
browser tab; they are not stored on the server or shared with other visitors. A deliberately fake example
to try is in [`dash_app/test_user_ensembles/test_user_ensemble.csv`](dash_app/test_user_ensembles/test_user_ensemble.csv)
(upload it as Greenland).

### Running the app locally

```bash
cd dash_app
pip install -r requirements.txt
python app.py              # then open http://localhost:8050
```

All the data the app needs are bundled in `dash_app/data/`, so no downloads happen at startup.

## The notebook

`exploring_fits/analyze_slr_predictions_interactive.ipynb` is the longer-form companion to the app,
and where most of the analysis was first developed. It:

- loads the ISMIP6 projections and IMBIE observations;
- reproduces and extends the figures of Aschwanden et al. (2021);
- includes interactive versions of the rate comparison;
- fits ANOVA and regression models of simulation misfit against model characteristics;
- compares the full-century Edwards et al. (2021) emulator projections with ISMIP6, checked against the
  paper's own published 2100 values.

The loaders it uses are in `utilities/`. They download the external datasets on first use and cache them
locally. The notebook needs a recent scientific Python stack: numpy, pandas, scipy, matplotlib, seaborn,
plotly, statsmodels, xarray, netCDF4, h5py, requests and ipywidgets.

## Data

| In the app | Source |
|---|---|
| Seroussi 2020 (ISMIP6 AIS) | Seroussi et al. (2020), *The Cryosphere* 14, 3033–3070, https://doi.org/10.5194/tc-14-3033-2020 |
| Goelzer 2020 (ISMIP6 GIS) | Goelzer et al. (2020), *The Cryosphere* 14, 3071–3096, https://doi.org/10.5194/tc-14-3071-2020 |
| DeConto & Pollard 2016 | DeConto & Pollard (2016), *Nature* 531, 591–597, https://doi.org/10.1038/nature17145; time series from Data Set S1 of Kopp et al. (2017), *Earth's Future* 5, 1217–1233, https://doi.org/10.1002/2017EF000663 (CC BY-NC-ND 4.0; bundled unmodified). Uncorrected runs, decadal values interpolated to annual |
| Aschwanden 2019 | Aschwanden et al. (2019), *Science Advances* 5, eaav9396, https://doi.org/10.1126/sciadv.aav9396; ensemble as archived with Aschwanden & Brinkerhoff (2022), https://doi.org/10.18739/A2KW57K4R |
| DeConto 2021 | DeConto et al. (2021), *Nature* 593, 83–89, https://doi.org/10.1038/s41586-021-03427-0; runs from the paper's Figure 1 source data (+1.5 °C, +2 °C, +3 °C and RCP8.5) |
| Coulon 2024 | Coulon et al. (2024), *The Cryosphere* 18, 653, https://doi.org/10.5194/tc-18-653-2024 |
| Rahlves 2025 | Rahlves et al. (2025), *The Cryosphere* 19, 1205, https://doi.org/10.5194/tc-19-1205-2025 |
| Goelzer 2025 (PROTECT GIS) | Goelzer et al. (2025), *The Cryosphere* 19, 6887, https://doi.org/10.5194/tc-19-6887-2025 |
| Edwards 2021 (AIS Main, AIS Risk Averse, GIS) | Edwards et al. (2021), *Nature* 593, 74–82, https://doi.org/10.1038/s41586-021-03302-y; samples from https://github.com/tamsinedwards/emulandice |
| IMBIE3 observations | Otosaka et al. (2026), *Scientific Data*, https://doi.org/10.1038/s41597-026-08088-0 |
| IMBIE2 observations | Otosaka et al. (2023), *Earth System Science Data* 15, 1597, https://doi.org/10.5194/essd-15-1597-2023 |

**All simulated mass changes are ice mass above flotation**, the part of an ice sheet that changes sea level
and that IMBIE's grounded-ice observations track. For ISMIP6 (both ice sheets) and Goelzer et al. (2025), mass
is ice volume above flotation (`ivaf`) × each model's own ice density (`rhoi`, 900–918 kg/m³). Rahlves et al.
(2025) provides mass above flotation directly. The studies published as sea-level contributions (Coulon,
Edwards, DeConto) are converted to mass at 362.5 Gt per mm, the same factor used by ISMIP6 and IMBIE. Aschwanden
et al. (2019) provides cumulative mass that its authors convert to sea level at the same factor.

`utilities/external_sources.py` documents how each external dataset was obtained and processed, and
the choices made along the way. If you use results from this tool, please cite the original studies
above.

## Repository layout

```
dash_app/        the web app: app.py (layout + callbacks), data.py, analysis.py, figures.py,
                 user_data.py, assets/, and the bundled data/
exploring_fits/  the notebook and the ISMIP6 files it reads
utilities/       data loaders shared by the notebook (and used to build dash_app/data/)
```

## Acknowledgments

This project began as, and was inspired by, Andy Aschwanden's
[ismip6-ipcc](https://github.com/aaschwanden/ismip6-ipcc) repository and its analysis of ISMIP6
projections for Aschwanden et al. (2021, *The Cryosphere*). Thanks also to the ISMIP6, IMBIE and
individual modeling teams whose openly archived results make this kind of comparison possible.

The app and much of the analysis code were developed with the help of Claude, Anthropic's AI assistant.
Tim Bartholomaus (University of Idaho) directed the work and takes responsibility for any errors in the
code or its results.

## Suggestions and corrections

Suggestions for new features, additional datasets, and corrections of any kind are very welcome. Please
[open an issue](https://github.com/tbartholomaus/ice_sheet_simulation_explorer/issues) or contact Tim
Bartholomaus.

## License

GNU General Public License v3.0. See [LICENSE](LICENSE).
