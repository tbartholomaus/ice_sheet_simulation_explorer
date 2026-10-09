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
from dash import Input, Output, State, ctx, dcc, html, no_update  # noqa: E402

import analysis as A  # noqa: E402
import figures as F  # noqa: E402
import user_data as U  # noqa: E402
from data import (  # noqa: E402
    ANOVA_CHARACTERISTICS, COMPOSITE_DIM, DEFAULT_OBS, EMULATION_SOURCES, ICE_SHEET_AREA_M2, DIM_LABEL, GROUP_DIMENSIONS, OBS_PRODUCTS, RUNS, SOURCE_COLOR, SOURCE_DEFAULT_CHECKED,
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

def _source_options():
    """Studies checklist with two sub-headings. A Checklist can't hold
    headings, so each is a disabled option (its checkbox hidden in
    style.css); its value never matches a study, so it can't affect a plot."""
    def study(s):
        return {"label": html.Span([html.Span(className="swatch", style={"background": SOURCE_COLOR[s]}), s]), "value": s}

    def heading(text, key):
        return {"label": html.Span(text, className="sources-heading"), "value": f"__heading_{key}", "disabled": True}

    physics = [s for s in SOURCE_LABELS if s not in EMULATION_SOURCES]
    emulation = [s for s in SOURCE_LABELS if s in EMULATION_SOURCES]
    return ([heading("Physics-based models", "physics")] + [study(s) for s in physics]
            + [heading("Emulations and syntheses", "emulation")] + [study(s) for s in emulation])


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


def graph(gid, height, **config):
    return dcc.Loading(type="dot", color="#0b63b6", children=dcc.Graph(
        id=gid, config={"displaylogo": False, "responsive": True, **config}, style={"height": f"{height}px"}))


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
    control("Simulation studies", html.Div([
        dcc.Checklist(id="sources", className="sources", value=list(SOURCE_DEFAULT_CHECKED), options=_source_options()),
        html.Div(id="user-chips"),  # one row per uploaded ensemble, with its own Remove link
        html.Button("+ Add your own ensemble", id="upload-open", className="upload-button", n_clicks=0),
    ])),
    control("Group by", html.Div([
        dcc.Dropdown(id="groupby", clearable=False, searchable=False, value="all", options=GROUP_OPTIONS),
        dcc.Checklist(id="collapse", className="toggle collapse-toggle", value=[],
                      options=[{"label": "Collapse RCPs and SSPs", "value": "on"}]),
        html.Div("RCP and SSP climate scenarios are not equivalent and can differ in their GMST change at 2100 "
                 "by as much as 0.5 °C.",
                 id="collapse-warning", className="warning", style={"display": "none"}),
    ]), "Collapsing combines simulations run under similar RCP and SSP scenarios (e.g. RCP8.5 and SSP5-8.5 "
        "become “composite 8.5”), setting aside differences in their climate trajectories."),
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
            "Drag to pan through time (out to 2100); scroll, use the box-zoom tool, or drag the y-axis to zoom, "
            "including vertically. Double-click to return to the window. "
            "Click a legend entry to hide it; double-click to show it alone.",
            [graph("ts-graph", 460, scrollZoom=True, doubleClick="autosize"),
             dcc.Store(id="ts-view"),  # the x range the user panned/zoomed to, for the current window
             html.H3("Mass change by 2100"),
             html.P("Change from 2015 to 2100 for each group (box: 25–75%; whiskers: 5–95%; line: median).",
                    className="lede"),
             # Compare projections from all simulations with only those whose
             # recent rate matches observations.
             html.Div(className="filter-bar", children=[
                 dcc.Checklist(id="match-filter", className="toggle", value=[],
                               options=[{"label": "Compare with only the simulations that match observations",
                                         "value": "on"}]),
                 html.Div(className="filter-tol", children=[
                     html.Span("Match: average rate over the averaging window within ±", className="filter-tol-label"),
                     html.Div(dcc.Slider(id="match-tol", min=5, max=100, step=5, value=25,
                                         marks={p: f"{p}%" for p in (5, 25, 50, 75, 100)},
                                         tooltip={"placement": "bottom", "template": "{value}%"}),
                              className="filter-tol-slider"),
                     html.Span("of the observed rate", className="filter-tol-label"),
                 ]),
                 html.Div(id="match-readout", className="filter-readout"),
             ]),
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
UPLOAD_EXAMPLE = """year,exp01,exp02,exp03
2015,0.0,0.0,0.0
2016,-231.4,-198.7,-260.2
2017,-455.0,-402.3,-512.8
..."""

upload_modal = html.Div(id="upload-modal", className="modal-backdrop", style={"display": "none"}, children=[
    html.Div(className="modal", role="dialog", children=[
        html.H2("Add your own ensemble"),
        html.P(["Upload a CSV file whose ", html.B("first column is the year"), " and whose other columns are "
                "your simulations, one per column, each holding ", html.B("cumulative mass (or mass change) in Gt"),
                " (mass loss negative). The header row names the columns: ", html.Code("year"),
                ", then an experiment name or number for each simulation."], className="lede"),
        html.Pre(UPLOAD_EXAMPLE, className="example"),
        html.P(f"Annual or finer time steps are fine. Up to {U.MAX_RUNS:,} simulations and "
               f"{U.MAX_BYTES // 1_000_000} MB per file, and up to {U.MAX_ENSEMBLES} ensembles at once, each "
               "shown as its own group. Your files stay in this browser tab; they aren't stored on the "
               "server or shared with anyone.", className="note"),
        dcc.Upload(id="upload-file", className="dropzone", accept=".csv,.txt,text/csv",
                   children=html.Div(["Drag a CSV here, or ", html.Span("choose a file", className="link")])),
        html.Div(id="upload-status", className="upload-status"),
        html.Div(id="upload-step2", style={"display": "none"}, children=[
            html.Div("Which ice sheet do these simulations represent?", className="control-label step-label"),
            dcc.RadioItems(id="upload-ice", className="segmented", inline=True, value=None,
                           options=[{"label": "Antarctica", "value": "AIS"}, {"label": "Greenland", "value": "GIS"}]),
            html.Div("Name to show in legends", className="control-label step-label"),
            dcc.Input(id="upload-name", type="text", value=U.DEFAULT_LABEL, maxLength=40, className="text-input"),
        ]),
        html.Div(className="modal-actions", children=[
            html.Button("Cancel", id="upload-cancel", className="button-secondary", n_clicks=0),
            html.Button("Add to plots", id="upload-add", className="button-primary", n_clicks=0, disabled=True),
        ]),
    ]),
])

app.layout = html.Div(className="page", children=[
    sidebar, main, upload_modal,
    dcc.Store(id="upload-pending"),                         # parsed, not yet confirmed
    dcc.Store(id="user-data", storage_type="session"),      # confirmed uploads, a list (this tab only)
])
server = app.server


# ── callbacks ─────────────────────────────────────────────────────────────

@app.callback(Output("years", "value"), Output("years-readout", "children"), Output("years", "max"),
              Input("years", "value"), Input("obs", "value"))
def _years(value, obs):
    if ctx.triggered_id == "obs":
        # Switching product: run the window up to the new product's last year
        # (2023 for IMBIE3), rather than leaving it where IMBIE2 capped it.
        value = [value[0], year_max(obs)]
    lo, hi = clamp_window(value, obs)
    fixed = [lo, hi] if [lo, hi] != list(value) else no_update
    return fixed, f"Jan {lo} – Dec {hi}", year_max(obs)


@app.callback(Output("rates-graph", "figure"),
              Input("years", "value"), Input("sources", "value"), Input("groupby", "value"),
              Input("units", "value"), Input("medians", "value"), Input("obs", "value"), Input("collapse", "value"),
              Input("user-data", "data"))
def _rates(years, sources, groupby, units, medians, obs, collapse, user):
    lo, hi = clamp_window(years, obs)
    return F.rates_figure(A.checked_mask(sources or []), lo, hi, _dim(groupby, bool(collapse)), units, bool(medians),
                          obs, U.from_store(user))


@app.callback(Output("ts-graph", "figure"),
              Input("years", "value"), Input("sources", "value"), Input("groupby", "value"), Input("units", "value"),
              Input("obs", "value"), Input("collapse", "value"), Input("user-data", "data"),
              State("ts-view", "data"))
def _timeseries(years, sources, groupby, units, obs, collapse, user, view):
    lo, hi = clamp_window(years, obs)
    valid = A.checked_mask(sources or [])
    # A rebuild for any control other than the window keeps the user's view;
    # an explicit range in a new figure would otherwise override uirevision.
    view_x = view["x"] if view and view.get("x") and (view.get("lo"), view.get("hi")) == (lo, hi) else None
    return F.timeseries_figure(valid, lo, hi, _dim(groupby, bool(collapse)), units, obs, U.from_store(user), view_x)


@app.callback(Output("y2100-graph", "figure"), Output("y2100-graph", "style"), Output("y2100-note", "children"),
              Output("match-readout", "children"),
              Input("years", "value"), Input("sources", "value"), Input("groupby", "value"), Input("units", "value"),
              Input("obs", "value"), Input("collapse", "value"), Input("user-data", "data"),
              Input("match-filter", "value"), Input("match-tol", "value"))
def _change2100(years, sources, groupby, units, obs, collapse, user, match_on, tol_pct):
    lo, hi = clamp_window(years, obs)
    valid = A.checked_mask(sources or [])
    gaps = []
    for pub in sorted(set(RUNS.loc[valid, "publication"])):
        m = valid & (RUNS["publication"] == pub).to_numpy()
        n_missing = int((~np.isfinite(A.CHANGE_2100[m])).sum())
        if n_missing:
            gaps.append(f"{n_missing} of {int(m.sum())} {pub} runs")
    users = U.from_store(user)
    for u in users:
        n_missing = int((~np.isfinite(U.change_2015_2100(u))).sum())
        if n_missing:
            gaps.append(f"{n_missing} of {len(u.vals)} {u.label} runs (need values in 2015 or 2016 and 2100)")
    note = f"Omitted (simulation ends before 2100): {'; '.join(gaps)}." if gaps else None
    dim = _dim(groupby, bool(collapse))
    matched, users_matched, readout = None, None, None
    if match_on:
        tol = (tol_pct or 25) / 100
        matched, bounds = A.obs_match(valid, lo, hi, obs, tol)
        users_matched = {}
        for u in users:
            r_obs, b_lo, b_hi = bounds[u.ice_sheet]
            r = U.rates(u, lo, hi)
            users_matched[u.label] = np.isfinite(r) & (r >= b_lo) & (r <= b_hi)
        readout = _match_readout(valid, matched, bounds, users, users_matched, units, obs, lo, hi, tol_pct)
    y2100 = F.change_2100_figure(valid, dim, units, users, matched, users_matched)
    return y2100, _height(y2100), note, readout


def _match_readout(valid, matched, bounds, users, users_matched, units, obs, lo, hi, tol_pct):
    """One line per ice sheet: the observed rate, the kept range, and how many
    simulations (with a 2015->2100 value) pass."""
    f, rate_label, _ = F.units_info(units)
    has2100 = np.isfinite(A.CHANGE_2100)
    lines = []
    for ice, name in (("AIS", "Antarctica"), ("GIS", "Greenland")):
        sel = valid & (_ICE == ice) & has2100
        n_all = int(sel.sum()) + sum(int(np.isfinite(U.change_2015_2100(u)).sum()) for u in users if u.ice_sheet == ice)
        if not n_all:
            continue
        n_kept = int((sel & matched).sum()) + sum(
            int((users_matched[u.label] & np.isfinite(U.change_2015_2100(u))).sum()) for u in users if u.ice_sheet == ice)
        r_obs, b_lo, b_hi = (v * f for v in bounds[ice])
        b_lo, b_hi = sorted((b_lo, b_hi))
        fmt = "{:.0f}" if units == "gt" else "{:.2f}"  # whole Gt/yr, or mm/yr to 0.01
        lines.append(html.Div([html.B(f"{name}: "), f"observed {fmt.format(r_obs)} {rate_label}; kept "
                               f"{fmt.format(b_lo)} to {fmt.format(b_hi)}: ", html.B(f"{n_kept} of {n_all}"),
                               " simulations."]))
    lines.append(html.Div(f"Rates over {lo}–{hi} against {OBS_PRODUCTS[obs]['label']}, within ±{tol_pct}%. "
                          "Faint outlined boxes: all selected simulations; solid boxes: those that match.",
                          className="muted"))
    return lines


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
              Input("collapse", "value"), Input("user-data", "data"))
def _bias(years, sources, groupby, units, scope, obs, collapse, user):
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
    uploads = []
    for u in U.from_store(user):
        if both or u.ice_sheet == scope:
            # Same definition as A.bias: rate minus the observed rate, per unit
            # area when pooling both ice sheets. Shown beside, never inside, the
            # ANOVA and best-match statistics -- an upload has no metadata.
            ub = U.rates(u, lo, hi) - A.imbie_rate(u.ice_sheet, lo, hi, obs)[0]
            uploads.append((u, ub * 1e12 / ICE_SHEET_AREA_M2[u.ice_sheet] if both else ub * F.units_info(units)[0]))
    if uploads:
        names = ", ".join(u.label for u, _ in uploads)
        verb = "is" if len(uploads) == 1 else "are"
        note = (note + " " if note else "") + (f"{names} {verb} shown for comparison only, not as part of the "
                                               "ANOVA or the best-match shares.")
    box = F.bias_figure(b, w, mask, dim, label, obs, uploads)
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


# ── upload your own ensemble ──────────────────────────────────────────────

@app.callback(Output("upload-modal", "style"),
              Input("upload-open", "n_clicks"), Input("upload-cancel", "n_clicks"), Input("upload-add", "n_clicks"),
              prevent_initial_call=True)
def _upload_modal(*_):
    return {"display": "flex" if ctx.triggered_id == "upload-open" else "none"}


@app.callback(Output("upload-pending", "data"), Output("upload-status", "children"),
              Output("upload-step2", "style"), Output("upload-status", "className"), Output("upload-name", "value"),
              Input("upload-file", "contents"), State("upload-file", "filename"), State("user-data", "data"),
              prevent_initial_call=True)
def _upload_parse(contents, filename, store):
    if not contents:
        return None, None, {"display": "none"}, "upload-status", no_update
    if len(U.entries(store)) >= U.MAX_ENSEMBLES:
        return (None, f"Up to {U.MAX_ENSEMBLES} ensembles can be shown at once; remove one first.",
                {"display": "none"}, "upload-status error", no_update)
    try:
        payload, summary = U.parse_upload(contents, filename)
    except ValueError as e:
        return None, f"{filename}: {e}", {"display": "none"}, "upload-status error", no_update
    # Default legend name: the file name, made unique among current uploads.
    stem = os.path.splitext(os.path.basename(filename or ""))[0][:40]
    name = U.unique_label(stem, {e["label"] for e in U.entries(store)})
    return payload, f"{filename}: {summary}", {"display": "block"}, "upload-status ok", name


@app.callback(Output("upload-add", "disabled"), Input("upload-pending", "data"), Input("upload-ice", "value"))
def _upload_ready(pending, ice):
    return not (pending and ice)


@app.callback(Output("user-data", "data"), Output("upload-file", "contents"), Output("upload-ice", "value"),
              Input("upload-add", "n_clicks"), Input({"type": "user-remove", "index": dash.ALL}, "n_clicks"),
              State("user-data", "data"), State("upload-pending", "data"), State("upload-ice", "value"),
              State("upload-name", "value"), prevent_initial_call=True)
def _upload_commit(_add, _removes, store, pending, ice, name):
    """Add the confirmed upload to the list, or drop the one whose Remove
    link was clicked. Each upload keeps the color it was given."""
    current = U.entries(store)
    if isinstance(ctx.triggered_id, dict):
        if not ctx.triggered[0]["value"]:  # chips re-rendered, not clicked
            return no_update, no_update, no_update
        label = ctx.triggered_id["index"]
        return [e for e in current if e["label"] != label], no_update, no_update
    if not (pending and ice) or len(current) >= U.MAX_ENSEMBLES:
        return no_update, no_update, no_update
    taken = {e["label"] for e in current}
    entry = {**pending, "ice_sheet": ice, "label": U.unique_label(name, taken),
             "color": U.free_color({e.get("color") for e in current})}
    # Reset the dialog so the next upload starts fresh.
    return current + [entry], None, None


@app.callback(Output("user-chips", "children"), Output("upload-open", "disabled"), Input("user-data", "data"))
def _user_chips(store):
    chips = []
    for u in U.from_store(store):
        n = len(u.names)
        where = {"AIS": "Antarctica", "GIS": "Greenland"}[u.ice_sheet]
        chips.append(html.Div(className="user-chip", children=[
            html.Span(className="swatch", style={"background": u.color}),
            html.Span(f"{u.label} ({where}, {n} simulation{'s' if n != 1 else ''})", className="user-chip-text"),
            html.Button("Remove", id={"type": "user-remove", "index": u.label},
                        className="link-button user-remove", n_clicks=0),
        ]))
    return chips, len(chips) >= U.MAX_ENSEMBLES


# After a sideways PAN of the time series (the visible span of years is
# unchanged), refit each panel's y axis to the data now in view: the series
# run to 2100, where values dwarf those near the window, and Plotly keeps the
# old y range while panning. Zooms are left alone -- a box zoom, scroll zoom
# or y-axis drag sets the vertical range the user asked for, and refitting
# would undo it. A double-click (which autoranges to 1950-2100) is turned
# into "back to the window". Runs in the browser; y-only relayouts are
# ignored, so it can't loop. The resulting x range goes to the ts-view store
# so a server rebuild keeps it.
app.clientside_callback(
    """
    function(relayout) {
        const nu = window.dash_clientside.no_update;
        if (!relayout || !Object.keys(relayout).some(k => k.startsWith("xaxis"))) return nu;
        const gd = document.querySelector("#ts-graph .js-plotly-plot");
        if (!gd || !gd._fullLayout || !gd._fullData) return nu;
        const meta = (gd.layout && gd.layout.meta) || {};
        const reset = relayout["xaxis.autorange"] || relayout["xaxis2.autorange"];
        const xr = reset && meta.window ? meta.window : gd._fullLayout.xaxis.range.slice();
        // Previous visible span, kept on window (the graph div gets replaced on
        // rebuilds) and keyed to the averaging window, which resets the view.
        const span = xr[1] - xr[0], key = meta.lo + "-" + meta.hi;
        const prev = window.__tsSpanKey === key && window.__tsSpan
            ? window.__tsSpan : (meta.window ? meta.window[1] - meta.window[0] : span);
        window.__tsSpan = span; window.__tsSpanKey = key;
        const panned = Math.abs(span - prev) < 1e-6 * Math.max(1, Math.abs(prev));
        const view = {lo: meta.lo, hi: meta.hi, x: reset ? null : xr};
        if (!reset && !panned) return view;  // a zoom: keep the user's vertical range
        const upd = {};
        if (reset && meta.window) upd["xaxis.range"] = meta.window.slice();
        [["x", "xaxis", "yaxis"], ["x2", "xaxis2", "yaxis2"]].forEach(([xid, xax, yax]) => {
            if (!gd._fullLayout[xax]) return;
            const [x0, x1] = xr;  // the panels share one x range (matches="x")
            let lo = Infinity, hi = -Infinity;
            gd._fullData.forEach(t => {
                if (t.xaxis !== xid || t.visible !== true || !t.x || !t.y) return;
                for (let i = 0; i < t.x.length; i++) {
                    const x = +t.x[i], y = +t.y[i];
                    if (x >= x0 && x <= x1 && isFinite(y)) { if (y < lo) lo = y; if (y > hi) hi = y; }
                }
            });
            if (isFinite(lo)) {
                const pad = 0.08 * Math.max(hi - lo, 1e-9);
                upd[yax + ".range"] = [lo - pad, hi + pad];
            }
        });
        if (Object.keys(upd).length) Plotly.relayout(gd, upd);
        return view;
    }
    """,
    Output("ts-view", "data"), Input("ts-graph", "relayoutData"), prevent_initial_call=True)


@app.callback(Output("collapse-warning", "style"), Input("collapse", "value"))
def _collapse_warning(collapse):
    return {"display": "block" if collapse else "none"}


if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8050)
