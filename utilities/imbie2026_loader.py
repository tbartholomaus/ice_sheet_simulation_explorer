import pandas as pd

from .helper import proj_start

# IMBIE 3 community assessment mass balance time series (Otosaka et al., 2026,
# Scientific Data, https://doi.org/10.1038/s41597-026-08088-0), covering
# 1971-2023 (Greenland) / 1979-2023 (Antarctica), from the NERC EDS UK Polar
# Data Centre (https://doi.org/10.5285/128c5e33-5224-4197-82f0-19dcc95b80a0,
# Open Government Licence v3.0, no restrictions).
#
# Successor to imbie2023_loader.py's Otosaka et al. (2023) ESSD product
# (that module is left as-is; this is a new module built from it, not a
# rewrite of it). The key difference driving that: each of these files
# already carries a self-consistent surface-mass-balance/dynamics
# partitioning (SMB + Dynamics == Total, confirmed directly against the raw
# data to ~1e-8 Gt/yr) for BOTH ice sheets, whereas imbie2023_loader.py's
# load_imbie2023_gis() had to merge in a separate, 2018-truncated
# Greenland-only dynamics workbook (imbie.org) that did NOT sum consistently
# with the main assessment on its own, and needed an unexplained
# "+/- 2*1964/10 Gt/yr" fudge factor to align the two. Antarctica had no
# partitioning at all under the old loader. Per explicit user decision
# (2026-09-22): use ONLY IMBIE3's own numbers (no old workbook, no fudge
# factor), for the two ice-sheet TOTALS only (not the WAIS/EAIS/Antarctic
# Peninsula regional breakdown IMBIE3 also publishes -- see
# IMBIE2026_GT_URLS's sibling files at the same DOI if that's wanted later).
IMBIE2026_GT_URLS = {
    "antarctica": (
        "https://ramadda.data.bas.ac.uk/repository/entry/get/imbie3_antarctica_Gt_partitioned.csv"
        "?entryid=synth:128c5e33-5224-4197-82f0-19dcc95b80a0:L2ltYmllM19hbnRhcmN0aWNhX0d0X3BhcnRpdGlvbmVkLmNzdg=="
    ),
    "greenland": (
        "https://ramadda.data.bas.ac.uk/repository/entry/get/imbie3_greenland_Gt_partitioned.csv"
        "?entryid=synth:128c5e33-5224-4197-82f0-19dcc95b80a0:L2ltYmllM19ncmVlbmxhbmRfR3RfcGFydGl0aW9uZWQuY3N2"
    ),
}


def _load_imbie2026(region):
    """
    Loads the IMBIE 3 community-assessment mass balance CSV for the given
    region ("antarctica" or "greenland") and normalizes it to this project's
    column names -- total mass change AND its SMB/dynamics partitioning,
    both sourced from this same file (no second source to merge in, unlike
    imbie2023_loader.py's GIS-only dynamics-workbook merge -- see this
    module's header comment).

    Each raw file opens with ~26 "# key: value" metadata lines (citation,
    license, time coverage, etc., all worth reading directly at the URLs in
    IMBIE2026_GT_URLS if there's ever a question about this data's
    provenance) before the real header row -- `comment="#"` skips those.
    Its own Date column is a monthly "YYYY-MM-01" string rather than a
    fractional Year; converted to a fractional-year Year column here to
    match every other loader in this project (misfit_dataframe,
    imbie_mass_loss_slope, plot_interactive_rate_comparison, etc. all key on
    a numeric Year) -- January of a given year lands on an exact `.0`, which
    is what the proj_start re-zeroing below relies on matching exactly.
    """
    df = pd.read_csv(IMBIE2026_GT_URLS[region], comment="#")
    date = pd.to_datetime(df["Date"])
    df["Year"] = date.dt.year + (date.dt.month - 1) / 12

    imbie = df.rename(
        columns={
            "Mass balance (Gt/yr)": "Rate of ice sheet mass change (Gt/yr)",
            "Cumulative mass balance anomaly (Gt)": "Cumulative ice sheet mass change (Gt)",
            "Cumulative mass balance anomaly uncertainty (Gt)": "Cumulative ice sheet mass change uncertainty (Gt)",
            "Surface mass balance anomaly (Gt/yr)": "Rate of surface mass balance anomaly (Gt/yr)",
            "Surface mass balance anomaly uncertainty (Gt/yr)": "Rate of surface mass balance anomaly uncertainty (Gt/yr)",
            "Cumulative surface mass balance anomaly (Gt)": "Cumulative surface mass balance anomaly (Gt)",
            "Cumulative surface mass balance anomaly uncertainty (Gt)": "Cumulative surface mass balance anomaly uncertainty (Gt)",
            "Dynamics mass balance anomaly (Gt/yr)": "Rate of ice dynamics anomaly (Gt/yr)",
            "Dynamics mass balance anomaly uncertainty (Gt/yr)": "Rate of ice dynamics anomaly uncertainty (Gt/yr)",
            "Cumulative dynamics mass balance anomaly (Gt)": "Cumulative ice dynamics anomaly (Gt)",
            "Cumulative dynamics mass balance anomaly uncertainty (Gt)": "Cumulative ice dynamics anomaly uncertainty (Gt)",
        }
    )[
        [
            "Year",
            "Cumulative ice sheet mass change (Gt)",
            "Cumulative ice sheet mass change uncertainty (Gt)",
            "Rate of ice sheet mass change (Gt/yr)",
            "Cumulative surface mass balance anomaly (Gt)",
            "Cumulative surface mass balance anomaly uncertainty (Gt)",
            "Rate of surface mass balance anomaly (Gt/yr)",
            "Rate of surface mass balance anomaly uncertainty (Gt/yr)",
            "Cumulative ice dynamics anomaly (Gt)",
            "Cumulative ice dynamics anomaly uncertainty (Gt)",
            "Rate of ice dynamics anomaly (Gt/yr)",
            "Rate of ice dynamics anomaly uncertainty (Gt/yr)",
        ]
    ].copy()

    # Re-zero every cumulative anomaly at proj_start (2015), same convention
    # imbie2023_loader.py uses -- everything downstream that plots a
    # cumulative series alongside ISMIP6 projections (which themselves start
    # at proj_start) expects a proj_start-relative-zero baseline.
    for v in [
        "Cumulative ice sheet mass change (Gt)",
        "Cumulative surface mass balance anomaly (Gt)",
        "Cumulative ice dynamics anomaly (Gt)",
    ]:
        imbie[v] -= imbie.loc[imbie["Year"] == proj_start, v].values

    # Same rebase-to-the-end-and-flip-sign transform imbie2023_loader.py
    # applies to the total's cumulative uncertainty (not applied to the
    # SMB/dynamics cumulative uncertainties there either, so not replicated
    # for them here). Preserved unchanged since nothing about switching data
    # sources should change this convention's meaning downstream.
    imbie["Cumulative ice sheet mass change uncertainty (Gt)"] -= imbie[
        "Cumulative ice sheet mass change uncertainty (Gt)"
    ].values[-1]
    imbie["Cumulative ice sheet mass change uncertainty (Gt)"] *= -1

    return imbie


def load_imbie2026_ais():
    """
    Loads the IMBIE 3 Antarctica mass balance time series (Otosaka et al.,
    2026, Scientific Data), including its native SMB/dynamics partitioning --
    unlike imbie2023_loader.py's load_imbie2023_ais(), which has no
    partitioning columns at all for Antarctica.
    """
    return _load_imbie2026("antarctica")


def load_imbie2026_gis():
    """
    Loads the IMBIE 3 Greenland mass balance time series (Otosaka et al.,
    2026, Scientific Data), including its native SMB/dynamics partitioning --
    unlike imbie2023_loader.py's load_imbie2023_gis(), no separate dynamics
    workbook is merged in here; see this module's header comment for why.
    """
    return _load_imbie2026("greenland")
