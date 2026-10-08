"""
Ice Sheet Simulation Explorer -- compares ice-sheet simulations (ISMIP6 and
other published ensembles) against IMBIE3 or IMBIE2 observations.

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
from dash import Input, Output, ctx, dcc, html, no_update  # noqa: E402

import analysis as A  # noqa: E402
import figures as F  # noqa: E402
from data import (  # noqa: E402
    ANOVA_CHARACTERISTICS, COMPOSITE_DIM, DEFAULT_OBS, DIM_LABEL, GROUP_DIMENSIONS, OBS_PRODUCTS, RUNS, SOURCE_COLOR, SOURCE_DEFAULT_CHECKED,
    SOURCE_LABELS,
)

YEAR_MIN, MIN_YEAR_SPAN = 2000, 5
YEAR_MAX = max(p["last_year"] for p in OBS_PRODUCTS.values())  # 2023, IMBIE3's last full year
YEAR_DEFAULT = [2015, 2023]
LOW_BIAS_FRAC = 0.10
_ICE = RUNS["ice_sheet"].to_numpy()


def year_max(obs):
    """Last year the window may reach: the chosen product's last full year."""
    return OBS_PRODUCTS.get(obs, OBS_PRODUCTS[DEFAULT_OBS])["last_year"]


def clamp_window(year_range, obs=DEFAULT_OBS):
    """Cut the window off at the observation product's last year, and widen a
    too-narrow window around its center, kept inside the slider range."""
    top = year_max(obs)
    lo, hi = year_range
    hi = min(hi, top)
    lo = min(lo, hi)
    if hi - lo >= MIN_YEAR_SPAN:
        return int(lo), int(hi)
    center = (lo + hi) / 2
    lo = max(YEAR_MIN, min(round(center - MIN_YEAR_SPAN / 2), top - MIN_YEAR_SPAN))
    return int(lo), int(lo + MIN_YEAR_SPAN)


GROUP_OPTIONS = [{"label": label, "value": key or "all"} for key, label in GROUP_DIMENSIONS]


def _dim(value, collapse=False):
    """Run-table column for a Group by value; with "Collapse RCPs and SSPs"
    on, the scenario grouping uses the pooled composite column instead."""
    if value in (None, "all"):
        return None
    return COMPOSITE_DIM if collapse and value == "scenario" else value


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
    control("Group by", html.Div([
        dcc.Dropdown(id="groupby", clearable=False, searchable=False, value="all", options=GROUP_OPTIONS),
        dcc.Checklist(id="collapse", className="toggle collapse-toggle", value=[],
                      options=[{"label": "Collapse RCPs and SSPs", "value": "on"}]),
        html.Div("RCP and SSP climate scenarios are not equivalent and can differ by as much as 0.5 °C.",
                 id="collapse-warning", className="warning", style={"display": "none"}),
    ]), "Collapsing pools each RCP with the SSP of the same 2100 forcing (e.g. RCP8.5 + SSP5-8.5 = "
        "“composite 8.5”) wherever climate scenarios are grouped or compared."),
    control("Units", dcc.RadioItems(
        id="units", className="segmented", value="gt", inline=True,
        options=[{"label": "Mass change", "value": "gt"}, {"label": "Sea level", "value": "sle"}])),
    control("Observations", dcc.RadioItems(
        id="obs", className="segmented", value=DEFAULT_OBS, inline=True,
        options=[{"label": p["label"], "value": k} for k, p in OBS_PRODUCTS.items()]),
        "IMBIE3 (2026) runs through 2023; IMBIE2 (2023) through 2020, so it caps the window at 2020."),
    control("Display", dcc.Checklist(
        id="medians", className="toggle", value=["on"], options=[{"label": "Show medians", "value": "on"}])),
    html.Div(className="sidebar-foot", children=[
        html.P(["Observations: ", *[x for k, p in OBS_PRODUCTS.items() for x in (
            f"{'; ' if k != next(iter(OBS_PRODUCTS)) else ''}{p['label']}, ",
            html.A(p["citation"], href=p["url"], target="_blank"))], "."]),
        html.P("Simulations: ISMIP6 (Seroussi et al. 2020; Goelzer et al. 2020) and other published ensembles. "
               "Large ensembles are weighted so each study counts about as much as one ISMIP6 institution."),
    ]),
])

main = html.Main(className="content", children=[
    html.Header(className="page-header", children=[
        html.H1("How do simulations of ice sheet mass loss compare with observations?"),
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
            "Median (line), 25–75% (darker band) and 5–95% (faint band) of simulated cumulative change, zeroed at "
            "the start of the averaging window (shaded), against IMBIE observations (black, hatched ±2σ). "
            "With more than 8 groups only the 25–75% band is drawn (hover for 5–95%). "
            "Click a legend entry to hide it; double-click to show it alone.",
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
             # A second Group by, kept in sync with the sidebar's, placed where
             # its effect on the plots below is visible.
             html.Div(className="group-picker", children=[
                 html.Label("Group by", htmlFor="bias-groupby", className="group-picker-label"),
                 dcc.Dropdown(id="bias-groupby", clearable=False, searchable=False, value="all",
                              options=GROUP_OPTIONS, className="group-picker-dropdown"),
                 html.Span("Choose which modeling characteristic to break the bias down by (or click a row "
                           "in the table). This is the same setting as Group by in the sidebar, so it also "
                           "regroups the plots above.", className="group-picker-hint"),
             ]),
             html.Div(className="two-col", children=[
                 html.Div([html.H3(id="bias-box-title"), html.Div(id="bias-box-note", className="note"),
                           graph("bias-graph", 430)]),
                 html.Div([html.H3(id="enrich-title"),
                           html.P("Share of each category among the best-matching simulations vs. among all "
                                  "selected simulations. Green bars longer than gray = over-represented "
                                  "among low-bias runs.", className="lede"),
                           graph("enrich-graph", 380)]),
             ])]),
    html.Footer(className="page-footer", children=[
        html.P("This web app was developed by Tim Bartholomaus at the University of Idaho, using claude ai."),
        html.P(["You can explore the source code for the app and an accompanying jupyter notebook at ",
                html.A("https://github.com/tbartholomaus/ice_sheet_simulation_explorer",
                       href="https://github.com/tbartholomaus/ice_sheet_simulation_explorer", target="_blank"),
                "."]),
    ]),
])

app = dash.Dash(__name__, title="Ice Sheet Simulation Explorer", external_stylesheets=[
    "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"])
app.layout = html.Div(className="page", children=[sidebar, main])
server = app.server


# ── callbacks ─────────────────────────────────────────────────────────────

@app.callback(Output("years", "value"), Output("years-readout", "children"), Output("years", "max"),
              Input("years", "value"), Input("obs", "value"))
def _years(value, obs):
    lo, hi = clamp_window(value, obs)
    fixed = [lo, hi] if [lo, hi] != list(value) else no_update
    return fixed, f"Jan {lo} – Dec {hi}", year_max(obs)


@app.callback(Output("rates-graph", "figure"),
              Input("years", "value"), Input("sources", "value"), Input("groupby", "value"),
              Input("units", "value"), Input("medians", "value"), Input("obs", "value"), Input("collapse", "value"))
def _rates(years, sources, groupby, units, medians, obs, collapse):
    lo, hi = clamp_window(years, obs)
    return F.rates_figure(A.checked_mask(sources or []), lo, hi, _dim(groupby, bool(collapse)), units, bool(medians), obs)


@app.callback(Output("ts-graph", "figure"), Output("y2100-graph", "figure"), Output("y2100-note", "children"),
              Input("years", "value"), Input("sources", "value"), Input("groupby", "value"), Input("units", "value"),
              Input("obs", "value"), Input("collapse", "value"))
def _timeseries(years, sources, groupby, units, obs, collapse):
    lo, hi = clamp_window(years, obs)
    valid = A.checked_mask(sources or [])
    gaps = []
    for pub in sorted(set(RUNS.loc[valid, "publication"])):
        m = valid & (RUNS["publication"] == pub).to_numpy()
        n_missing = int((~np.isfinite(A.CHANGE_2100[m])).sum())
        if n_missing:
            gaps.append(f"{n_missing} of {int(m.sum())} {pub} runs")
    note = f"Omitted (simulation ends before 2100): {'; '.join(gaps)}." if gaps else None
    dim = _dim(groupby, bool(collapse))
    return (F.timeseries_figure(valid, lo, hi, dim, units, obs), F.change_2100_figure(valid, dim, units), note)


def _fmt_p(p):
    if not np.isfinite(p):
        return "—"
    return "< 0.001" if p < 0.001 else f"{p:.3f}"


def _anova_table(one_way, combined, shown=None):
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
        c = r["characteristic"]
        name = [html.Button(DIM_LABEL[c], id={"type": "anova-pick", "dim": c}, className="link-button",
                            title=f"Show bias by {_lower(DIM_LABEL[c])} below"),
                html.Span(f" ({r['n_categories']})", className="muted")]
        if c == shown:
            name.append(html.Span("shown below", className="tag"))
        rows.append(html.Tr(className="current" if c == shown else None, children=[
            html.Td(name),
            html.Td(html.Div(className="bar-cell", children=[
                html.Div(className="bar", style={"width": f"{100 * r['r2']:.1f}%"}),
                html.Span(f"{100 * r['r2']:.0f}%")])),
            html.Td(_fmt_p(r["p"]), className="sig" if r["p"] < 0.05 else ""),
            *joint]))
    return html.Table(className="anova", children=[html.Thead(head), html.Tbody(rows)])


def _lower(label):
    return " ".join(w if w.isupper() else w.lower() for w in label.split())  # keep "GCM", "RCP/SSP"


def _height(fig, default=300):
    return {"height": f"{(fig.layout.height if fig is not None else None) or default}px"}


@app.callback(Output("bias-summary", "children"), Output("anova-table", "children"),
              Output("bias-graph", "figure"), Output("bias-box-title", "children"), Output("bias-box-note", "children"),
              Output("enrich-graph", "figure"), Output("enrich-title", "children"),
              # Both figures grow with their number of categories; the
              # container has to grow with them or the rows get squashed.
              Output("bias-graph", "style"), Output("enrich-graph", "style"),
              Input("years", "value"), Input("sources", "value"), Input("groupby", "value"),
              Input("units", "value"), Input("bias-scope", "value"), Input("obs", "value"),
              Input("collapse", "value"))
def _bias(years, sources, groupby, units, scope, obs, collapse):
    lo, hi = clamp_window(years, obs)
    valid = A.checked_mask(sources or []) & np.isfinite(A.run_rates(lo, hi))
    both = scope == "both"
    mask = valid & (np.isin(_ICE, ["AIS", "GIS"]) if both else _ICE == scope)
    chars = [COMPOSITE_DIM if collapse and c == "scenario" else c for c in ANOVA_CHARACTERISTICS]
    if both:
        # Only categories present for BOTH ice sheets can be compared; anything
        # else would just re-encode which ice sheet a run belongs to.
        shared = {c: set(RUNS.loc[mask & (_ICE == "AIS"), c]) & set(RUNS.loc[mask & (_ICE == "GIS"), c]) for c in chars}
        chars = [c for c in chars if len(shared[c]) >= 2]
        for c in chars:
            mask &= RUNS[c].isin(shared[c]).to_numpy()
        label = "mm/yr water equivalent"
        b = A.bias(mask, lo, hi, obs, area_normalized=True)
    else:
        f, label, _ = F.units_info(units)
        b = A.bias(mask, lo, hi, obs) * f
    w = A.run_weights(mask)
    empty = F.go.Figure(layout=dict(template=F.TEMPLATE))
    if mask.sum() < 3:
        msg = html.P("Not enough simulations selected for this comparison.", className="muted")
        return msg, None, empty, "", "", empty, "", _height(None), _height(None)

    one_way_chars = chars + (["publication"] if len(set(RUNS.loc[mask, "publication"])) >= 2 else [])
    factors = {c: RUNS[c].to_numpy()[mask] for c in one_way_chars}
    one_way, _ = A.anova(b[mask], w[mask], factors)
    _, combined = A.anova(b[mask], w[mask], {c: factors[c] for c in chars})
    if not one_way:
        return html.P("No characteristic varies across the selected simulations.", className="muted"), None, \
            empty, "", "", empty, "", _height(None), _height(None)

    top = one_way[0]
    summary = [html.Span("Strongest single factor: "), html.B(DIM_LABEL[top["characteristic"]]),
               f" explains {100 * top['r2']:.0f}% of the spread in bias (p {_fmt_p(top['p'])})."]
    if np.isfinite(combined["r2"]):
        summary.append(f" Together the modeling characteristics explain {100 * combined['r2']:.0f}%.")
    summary.append(f" Based on {int(mask.sum())} simulations, {lo}–{hi}, against {OBS_PRODUCTS[obs]['label']}.")

    dim, note = _dim(groupby, bool(collapse)), ""
    if dim is None or dim not in one_way_chars:
        if dim is None:
            note = "“All simulations” has no categories to compare, so the strongest factor is shown."
        else:
            note = (f"{DIM_LABEL[dim]} doesn't vary across this selection (or isn't shared by both ice "
                    "sheets), so the strongest factor is shown.")
        dim = top["characteristic"]
    box = F.bias_figure(b, w, mask, dim, label, obs)
    rows, cut = A.low_bias_enrichment(b[mask], w[mask], RUNS[dim].to_numpy()[mask], LOW_BIAS_FRAC)
    enrich = F.enrichment_figure(rows, dim, LOW_BIAS_FRAC)
    return (summary, _anova_table(one_way, combined, dim), box, f"Bias by {_lower(DIM_LABEL[dim])}", note, enrich,
            f"What do the best-matching runs share? (|bias| ≤ {cut:.3g} {label})", _height(box), _height(enrich))


@app.callback(Output("groupby", "value"), Output("bias-groupby", "value"),
              Input("groupby", "value"), Input("bias-groupby", "value"),
              Input({"type": "anova-pick", "dim": dash.ALL}, "n_clicks"), prevent_initial_call=True)
def _sync_groupby(sidebar, local, _clicks):
    """One Group by setting, three ways to change it: the sidebar dropdown,
    the dropdown in the bias section, or a click on an ANOVA table row."""
    src = ctx.triggered_id
    if isinstance(src, dict):
        if not ctx.triggered[0]["value"]:  # table re-rendered, not clicked
            return no_update, no_update
        dim = "scenario" if src["dim"] == COMPOSITE_DIM else src["dim"]
        return dim, dim
    if src == "bias-groupby":
        return local, no_update
    return no_update, sidebar


@app.callback(Output("collapse-warning", "style"), Input("collapse", "value"))
def _collapse_warning(collapse):
    return {"display": "block" if collapse else "none"}


if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8050)
