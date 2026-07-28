import pandas as pd

from .helper import proj_start

# IMBIE community assessment mass balance time series (Otosaka et al., 2023, ESSD,
# https://doi.org/10.5194/essd-15-1597-2023), 1992-2020, from the BAS Polar Data Centre
# (https://doi.org/10.5285/77B64C55-7166-4A06-9DEF-2E400398E452).
IMBIE2023_GT_URLS = {
    "antarctica": (
        "https://ramadda.data.bas.ac.uk/repository/entry/get/imbie_antarctica_2021_Gt.csv"
        "?entryid=synth:77b64c55-7166-4a06-9def-2e400398e452:L2ltYmllX2FudGFyY3RpY2FfMjAyMV9HdC5jc3Y="
    ),
    "greenland": (
        "https://ramadda.data.bas.ac.uk/repository/entry/get/imbie_greenland_2021_Gt.csv"
        "?entryid=synth:77b64c55-7166-4a06-9def-2e400398e452:L2ltYmllX2dyZWVubGFuZF8yMDIxX0d0LmNzdg=="
    ),
}


def _load_imbie2023(region):
    """
    Loads the IMBIE 2023 community-assessment mass balance CSV for the given
    region ("antarctica" or "greenland") and normalizes it to this notebook's
    column names.
    """
    df = pd.read_csv(IMBIE2023_GT_URLS[region])
    imbie = df.rename(
        columns={
            "Mass balance (Gt/yr)": "Rate of ice sheet mass change (Gt/yr)",
            "Cumulative mass balance (Gt)": "Cumulative ice sheet mass change (Gt)",
            "Cumulative mass balance uncertainty (Gt)": "Cumulative ice sheet mass change uncertainty (Gt)",
        }
    )[
        [
            "Year",
            "Cumulative ice sheet mass change (Gt)",
            "Cumulative ice sheet mass change uncertainty (Gt)",
            "Rate of ice sheet mass change (Gt/yr)",
        ]
    ].copy()

    imbie["Cumulative ice sheet mass change (Gt)"] -= imbie.loc[
        imbie["Year"] == proj_start, "Cumulative ice sheet mass change (Gt)"
    ].values

    imbie["Cumulative ice sheet mass change uncertainty (Gt)"] -= imbie[
        "Cumulative ice sheet mass change uncertainty (Gt)"
    ].values[-1]
    imbie["Cumulative ice sheet mass change uncertainty (Gt)"] *= -1

    return imbie


def load_imbie2023_ais():
    """
    Loading the IMBIE 2023 Antarctica mass balance time series
    (Otosaka et al., 2023, ESSD).
    """
    return _load_imbie2023("antarctica")


def load_imbie2023_gis():
    """
    Loading the IMBIE 2023 Greenland mass balance time series
    (Otosaka et al., 2023, ESSD), merged with the SMB/ice-dynamics
    partitioning from the original Greenland-specific dynamics workbook
    (http://imbie.org/wp-content/uploads/2012/11/imbie_dataset_greenland_dynamics-2020_02_28.xlsx),
    since the 2023 assessment product does not include that decomposition.
    The dynamics workbook only extends to 2018.9, so partitioning columns
    are NaN beyond that even though the headline mass-balance columns
    extend to 2020.9.
    """
    imbie = _load_imbie2023("greenland")

    dyn_df = pd.read_excel(
        "http://imbie.org/wp-content/uploads/2012/11/imbie_dataset_greenland_dynamics-2020_02_28.xlsx",
        sheet_name="Greenland Ice Mass",
        engine="openpyxl",
    )[
        [
            "Year",
            "Cumulative surface mass balance anomaly (Gt)",
            "Cumulative surface mass balance anomaly uncertainty (Gt)",
            "Cumulative ice dynamics anomaly (Gt)",
            "Cumulative ice dynamics anomaly uncertainty (Gt)",
            "Rate of mass balance anomaly (Gt/yr)",
            "Rate of ice dynamics anomaly (Gt/yr)",
            "Rate of mass balance anomaly uncertainty (Gt/yr)",
            "Rate of ice dyanamics anomaly uncertainty (Gt/yr)",
        ]
    ].rename(
        columns={
            "Rate of mass balance anomaly (Gt/yr)": "Rate of surface mass balance anomaly (Gt/yr)",
            "Rate of mass balance anomaly uncertainty (Gt/yr)": "Rate of surface mass balance anomaly uncertainty (Gt/yr)",
            "Rate of ice dyanamics anomaly uncertainty (Gt/yr)": "Rate of ice dynamics anomaly uncertainty (Gt/yr)",
        }
    ).copy()

    for v in [
        "Cumulative ice dynamics anomaly (Gt)",
        "Cumulative surface mass balance anomaly (Gt)",
    ]:
        dyn_df[v] -= dyn_df.loc[dyn_df["Year"] == proj_start, v].values

    dyn_df["Rate of surface mass balance anomaly (Gt/yr)"] += 2 * 1964 / 10
    dyn_df["Rate of ice dynamics anomaly (Gt/yr)"] -= 2 * 1964 / 10

    # Merge on rounded Year: both sources are monthly grids (steps of 1/12) but the
    # dynamics workbook stores more decimal digits, so exact float equality can miss.
    imbie["_yr_round"] = imbie["Year"].round(4)
    dyn_df["_yr_round"] = dyn_df["Year"].round(4)
    imbie = pd.merge(imbie, dyn_df.drop(columns="Year"), on="_yr_round", how="left").drop(columns="_yr_round")

    return imbie
