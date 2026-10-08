"""
Ice Sheet Simulation Explorer -- compares ice-sheet simulations (ISMIP6 and
other published ensembles) against IMBIE3 observations.

Layout: a sidebar of shared controls and one scrolling page of sections --
rates, mass change through time, mass change by 2100, and what drives bias
(weighted ANOVA). Every control is a plain Dash input; each section has its
own callback so it updates (and shows its spinner) independently.

Modules: data.py (loading + run table), analysis.py (numerics),
figures.py (plots), this file (layout + callbacks).

Run locally:   pip install -r requirements.txt && python app.py
Deploy:        gunicorn app:server --timeout 120   (see render.yaml)
"""

import os

# Single-threaded BLAS/OpenMP, set before numpy is imported -- on Render's
# thin CPU quota, oversized thread pools caused a ~100x slowdown.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import dash  # noqa: E402
from dash import Input, Output, dcc, html, no_update  # noqa: E402

import analysis as A  # noqa: E402
import figures as F  # noqa: E402
from data import (  # noqa: E402
    ANOVA_CHARACTERISTICS, DIM_LABEL, GROUP_DIMENSIONS, RUNS, SOURCE_COLOR, SOURCE_DEFAULT_CHECKED,
    SOURCE_LABELS,
)

YEAR_MIN, YEAR_MAX, MIN_YEAR_SPAN = 2000, 2023, 5  # 2023 == IMBIE3's last full year
YEAR_DEFAULT = [2015, 2023]
LOW_BIAS_FRAC = 0.10
_ICE = RUNS["ice_sheet"].to_numpy()


def clamp_window(year_range):
    """Widen a too-narrow window around its center, kept inside the slider range."""
    lo, hi = year_range
    if hi - lo >= MIN_YEAR_SPAN:
        return int(lo), int(hi)
    center = (lo + hi) / 2
    lo = max(YEAR_MIN, min(round(center - MIN_YEAR_SPAN / 2), YEAR_MAX - MIN_YEAR_SPAN))
    return int(lo), int(lo + MIN_YEAR_SPAN)


def _dim(value):
    return None if value in (None, "all") else value


# ── layout ────────────────────────────────────────────────────────────────

def control(label, child, hint=None):
    return html.Div(className="control", children=[
        html.Div(label, className="control-label"),
        child,
        html.Div(hint, className="control-hint") if hint else None,
    ])


def section(sid, title, lede, children):
    return html.Section(id=sid, className="card", children=[
        html.H2(title), html.P(lede, className="lede"), *children,
    ])


def graph(gid, height):
    return dcc.Loading(type="dot", color="#0b63b6", children=dcc.Graph(
        id=gid, config={"displaylogo": False, "responsive": True}, style={"height": f"{height}px"}))


sidebar = html.Aside(className="sidebar", children=[
    html.Div(className="brand", children=[
        html.Div("Ice Sheet", className="brand-top"),
        html.Div("Simulation Explorer", className="brand-bottom"),
    ]),
    control("Averaging window", html.Div([
        dcc.RangeSlider(id="years", min=YEAR_MIN, max=YEAR_MAX, step=1, value=YEAR_DEFAULT, allowCross=False,
                        marks={y: str(y) for y in (2000, 2005, 2010, 2015, 2020)},
                        tooltip={"placement": "bottom"}),
        html.Div(id="years-readout", className="readout"),
    ]), "Rates and bias are averaged over this window; time series are zeroed at its start."),
    control("Simulation studies", dcc.Checklist(
        id="sources", className="sources", value=list(SOURCE_DEFAULT_CHECKED),
        options=[{"label": html.Span([html.Span(className="swatch", style={"background": SOURCE_COLOR[s]}), s]),
                  "value": s} for s in SOURCE_LABELS])),
    control("Group by", dcc.Dropdown(
        id="groupby", clearable=False, searchable=False, value="all",
        options=[{"label": label, "value": key or "all"} for key, label in GROUP_DIMENSIONS])),
    control("Units", dcc.RadioItems(
        id="units", className="segmented", value="gt", inline=True,
        options=[{"label": "Mass change", "value": "gt"}, {"label": "Sea level", "value": "sle"}])),
    control("Display", dcc.Checklist(
        id="medians", className="toggle", value=["on"], options=[{"label": "Show medians", "value": "on"}])),
    html.Div(className="sidebar-foot", children=[
        html.P(["Observations: IMBIE3, ", html.A("Otosaka et al. (2026)", href="https://doi.org/10.1038/s41597-026-08088-0",
                                                 target="_blank"), "."]),
        html.P("Simulations: ISMIP6 (Seroussi et al. 2020; Goelzer et al. 2020) and other published ensembles. "
               "Large ensembles are weighted so each study counts about as much as one ISMIP6 institution."),
    ]),
])

main = html.Main(className="content", children=[
    html.Header(className="page-header", children=[
        html.H1("How well do ice sheet simulations match observed mass loss?"),
        html.P("Compare simulated rates of ice sheet mass change with satellite observations, see how "
               "simulations evolve through time and by 2100, and explore which modeling choices drive "
               "the differences. Hover for details; drag to zoom; double-click to reset.",
               className="lede"),
    ]),
    section("rates", "Rates of mass change",
            "Distribution of each simulation's average rate over the window, against the observed IMBIE rate. "
            "Each dot is one simulation; curves are weighted density estimates; triangles mark medians.",
            [graph("rates-graph", 760)]),
    section("timeseries", "Mass change through time",
            "Median (line) and 5–95% range (band) of simulated cumulative change, zeroed at the start of the "
            "averaging window (shaded), against IMBIE observations.",
            [graph("ts-graph", 460),
             html.H3("Mass change by 2100"),
             html.P("Change from 2015 to 2100 for each group (box: 25–75%; whiskers: 5–95%; line: median).",
                    className="lede"),
             graph("y2100-graph", 440),
             html.Div(id="y2100-note", className="note")]),
    section("bias", "What drives bias?",
            "Bias is each simulation's rate minus the observed rate over the averaging window. The table shows "
            "how much of the spread in bias each modeling characteristic explains (weighted ANOVA).",
            [html.Div(className="bias-toolbar", children=[
                dcc.RadioItems(id="bias-scope", className="segmented", value="AIS", inline=True,
                               options=[{"label": "Antarctica", "value": "AIS"},
                                        {"label": "Greenland", "value": "GIS"},
                                        {"label": "Both (per unit area)", "value": "both"}]),
             ]),
             html.Div(id="bias-summary", className="summary"),
             dcc.Loading(type="dot", color="#0b63b6", children=html.Div(id="anova-table")),
             html.Div(className="two-col", children=[
                 html.Div([html.H3(id="bias-box-title"), graph("bias-graph", 430)]),
                 html.Div([html.H3(id="enrich-title"),
                           html.P("Share of each category among the best-matching simulations vs. among all "
                                  "selected simulations. Green bars longer than gray = over-represented "
                                  "among low-bias runs.", className="lede"),
                           graph("enrich-graph", 380)]),
             ])]),
])

app = dash.Dash(__name__, title="Ice Sheet Simulation Explorer", external_stylesheets=[
    "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"])
app.layout = html.Div(className="page", children=[sidebar, main])
server = app.server


# ── callbacks ─────────────────────────────────────────────────────────────

@app.callback(Output("years", "value"), Output("years-readout", "children"), Input("years", "value"))
def _years(value):
    lo, hi = clamp_window(value)
    fixed = [lo, hi] if [lo, hi] != list(value) else no_update
    return fixed, f"Jan {lo} – Dec {hi}"


@app.callback(Output("rates-graph", "figure"),
              Input("years", "value"), Input("sources", "value"), Input("groupby", "value"),
              Input("units", "value"), Input("medians", "value"))
def _rates(years, sources, groupby, units, medians):
    lo, hi = clamp_window(years)
    return F.rates_figure(A.checked_mask(sources or []), lo, hi, _dim(groupby), units, bool(medians))


@app.callback(Output("ts-graph", "figure"), Output("y2100-graph", "figure"), Output("y2100-note", "children"),
              Input("years", "value"), Input("sources", "value"), Input("groupby", "value"), Input("units", "value"))
def _timeseries(years, sources, groupby, units):
    lo, hi = clamp_window(years)
    valid = A.checked_mask(sources or [])
    gaps = []
    for pub in sorted(set(RUNS.loc[valid, "publication"])):
        m = valid & (RUNS["publication"] == pub).to_numpy()
        n_missing = int((~np.isfinite(A.CHANGE_2100[m])).sum())
        if n_missing:
            gaps.append(f"{n_missing} of {int(m.sum())} {pub} runs")
    note = f"Omitted (simulation ends before 2100): {'; '.join(gaps)}." if gaps else None
    return (F.timeseries_figure(valid, lo, hi, _dim(groupby), units),
            F.change_2100_figure(valid, _dim(groupby), units), note)


def _fmt_p(p):
    if not np.isfinite(p):
        return "—"
    return "< 0.001" if p < 0.001 else f"{p:.3f}"


def _anova_table(one_way, combined):
    terms = {t["characteristic"]: t for t in combined["terms"]}
    head = html.Tr([html.Th("Characteristic"), html.Th("Variance explained alone"), html.Th("p"),
                    html.Th("Unique effect with all others (F)"), html.Th("p")])
    rows = []
    for r in one_way:
        t = terms.get(r["characteristic"])
        if t is None:
            joint = [html.Td("not in joint model", className="muted", colSpan=2)]
        elif t["df"] == 0:
            joint = [html.Td("confounded with others", className="muted", colSpan=2)]
        else:
            joint = [html.Td(f"{t['F']:.1f}"), html.Td(_fmt_p(t["p"]), className="sig" if t["p"] < 0.05 else "")]
        rows.append(html.Tr([
            html.Td(f"{DIM_LABEL[r['characteristic']]} ({r['n_categories']})"),
            html.Td(html.Div(className="bar-cell", children=[
                html.Div(className="bar", style={"width": f"{100 * r['r2']:.1f}%"}),
                html.Span(f"{100 * r['r2']:.0f}%")])),
            html.Td(_fmt_p(r["p"]), className="sig" if r["p"] < 0.05 else ""),
            *joint]))
    return html.Table(className="anova", children=[html.Thead(head), html.Tbody(rows)])


@app.callback(Output("bias-summary", "children"), Output("anova-table", "children"),
              Output("bias-graph", "figure"), Output("bias-box-title", "children"),
              Output("enrich-graph", "figure"), Output("enrich-title", "children"),
              Input("years", "value"), Input("sources", "value"), Input("groupby", "value"),
              Input("units", "value"), Input("bias-scope", "value"))
def _bias(years, sources, groupby, units, scope):
    lo, hi = clamp_window(years)
    valid = A.checked_mask(sources or []) & np.isfinite(A.run_rates(lo, hi))
    both = scope == "both"
    mask = valid & (np.isin(_ICE, ["AIS", "GIS"]) if both else _ICE == scope)
    chars = list(ANOVA_CHARACTERISTICS)
    if both:
        # Only categories present for BOTH ice sheets can be compared; anything
        # else would just re-encode which ice sheet a run belongs to.
        shared = {c: set(RUNS.loc[mask & (_ICE == "AIS"), c]) & set(RUNS.loc[mask & (_ICE == "GIS"), c]) for c in chars}
        chars = [c for c in chars if len(shared[c]) >= 2]
        for c in chars:
            mask &= RUNS[c].isin(shared[c]).to_numpy()
        label = "mm/yr water equivalent"
        b = A.bias(mask, lo, hi, area_normalized=True)
    else:
        f, label, _ = F.units_info(units)
        b = A.bias(mask, lo, hi) * f
    w = A.run_weights(mask)
    empty = F.go.Figure(layout=dict(template=F.TEMPLATE))
    if mask.sum() < 3:
        msg = html.P("Not enough simulations selected for this comparison.", className="muted")
        return msg, None, empty, "", empty, ""

    one_way_chars = chars + (["publication"] if len(set(RUNS.loc[mask, "publication"])) >= 2 else [])
    factors = {c: RUNS[c].to_numpy()[mask] for c in one_way_chars}
    one_way, _ = A.anova(b[mask], w[mask], factors)
    _, combined = A.anova(b[mask], w[mask], {c: factors[c] for c in chars})
    if not one_way:
        return html.P("No characteristic varies across the selected simulations.", className="muted"), None, \
            empty, "", empty, ""

    top = one_way[0]
    summary = [html.Span("Strongest single factor: "), html.B(DIM_LABEL[top["characteristic"]]),
               f" explains {100 * top['r2']:.0f}% of the spread in bias (p {_fmt_p(top['p'])})."]
    if np.isfinite(combined["r2"]):
        summary.append(f" Together the modeling characteristics explain {100 * combined['r2']:.0f}%.")
    summary.append(f" Based on {int(mask.sum())} simulations, {lo}–{hi}.")

    dim = _dim(groupby)
    if dim is None or dim not in one_way_chars:
        dim = top["characteristic"]
    box = F.bias_figure(b, w, mask, dim, label)
    rows, cut = A.low_bias_enrichment(b[mask], w[mask], RUNS[dim].to_numpy()[mask], LOW_BIAS_FRAC)
    enrich = F.enrichment_figure(rows, dim, LOW_BIAS_FRAC)
    return (summary, _anova_table(one_way, combined), box, f"Bias by {DIM_LABEL[dim].lower()}", enrich,
            f"What do the best-matching runs share? (|bias| ≤ {cut:.3g} {label})")


if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8050)
