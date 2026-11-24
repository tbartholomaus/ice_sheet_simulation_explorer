# Simulation metadata

Everything the app and the notebook know about each simulation — which ice sheet model, initialization,
sliding law, stress balance, resolution, climate model, scenario and ocean forcing it used — is in these
CSV files. Nothing is hard-coded elsewhere, so the encoding can be checked here without reading any code.

**Please check these encodings**, especially for models or studies you know well. If you disagree with one,
[open an issue](https://github.com/tbartholomaus/ice_sheet_simulation_explorer/issues) naming the file, the
row and your source.
Lines starting with `#` are comments (citations and notes); spreadsheet programs may need to skip them.

## Files

| File | What it holds |
|---|---|
| `<Study>_models.csv` | One row per model configuration (Group, Model). Every configuration is independent, including different versions of the same model. |
| `<Study>_experiments.csv` | One row per experiment. Where the study has an experiment table, it is transcribed verbatim (the paper's own column names and values). |
| `experiment_columns.csv` | For each study, which experiments column fills which app field, and how verbatim values translate to encoded ones (e.g. `Medium*` → `Medium`). |
| `encoding_rules.csv` | Every allowed encoded value of each field, in display order, with its meaning and the rule that assigns it. |

Studies: `Seroussi2020` (ISMIP6 Antarctica), `Goelzer2020` (ISMIP6 Greenland), `DeConto2016`,
`Aschwanden2019`, `DeConto2021`, `Coulon2024`, `Rahlves2025`, `Goelzer2025`, `Edwards2021`, `IPCCAR6`.

## Reading a models file

For each characteristic there is a group of columns:

| Column | Meaning |
|---|---|
| `<field>_described` | How the source describes it. Text in double quotes is quoted from the source. |
| `<field>_source` | Where that description comes from (paper, table, appendix section). |
| `<field>_encoded` | The value the app and notebook use for grouping and the ANOVA. |
| `<field>_note` | Ambiguities, disagreements between sources, and how they were resolved. |

`<field>` is `ice_model`, `stress_balance`, `resolution`, `initialization` or `sliding_law`. Supporting
columns: `stress_balance_code` and `initialization_code` (the source's own codes), `res_min_km` /
`res_max_km` (the resolution bin is set from `res_min_km`), `initial_year`, `sliding_exponent`,
`effective_pressure`, `sliding_law_confidence`, and the open/standard melt parameterizations.

`review_status` is `reviewed`, `not yet reviewed`, or `not applicable` (emulators and assessments, which
have no ice sheet model physics). For studies not yet reviewed, the `*_encoded` values are the encodings in
use today, marked `(unreviewed)` or `Not yet reviewed` so they never pool with reviewed values; their
`*_described` columns are filled in when the study is reviewed.

## Review reports

Each reviewed study has a report in [`docs/metadata_review/`](../../docs/metadata_review/) giving the sources
consulted, every discrepancy found, and the reasoning and date of each encoding decision:
[`01_Seroussi2020_ISMIP6_AIS.md`](../../docs/metadata_review/01_Seroussi2020_ISMIP6_AIS.md) and
[`02_Goelzer2020_ISMIP6_GIS.md`](../../docs/metadata_review/02_Goelzer2020_ISMIP6_GIS.md).

## Checks

`dash_app/metadata.py` validates every file when the app or notebook starts: each encoded value must be
allowed by `encoding_rules.csv`, each resolution bin must match `res_min_km`, every experiments column must
be described in `experiment_columns.csv`, and keys must be unique. `dash_app/tests/test_metadata.py` also
checks that every bundled simulation has a models row and an experiments row.
