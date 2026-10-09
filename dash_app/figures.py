"""Figure builders -- one per page section, each carrying only the active
grouping (unlike the original single figure that pre-built every grouping
for client-side toggling, which is what drove the 67 MB responses)."""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import analysis as A
import user_data as U
from data import COLOR_MAPS, DEFAULT_OBS, DIM_LABEL, OBS_PRODUCTS, GT_TO_MM_SLE, ICE_SHEET_NAMES, ICE_SHEETS, RUNS

FONT = '"Inter", "Segoe UI", system-ui, -apple-system, sans-serif'
INK = "#1f2a3c"
MUTED = "#6b7585"
GRID = "#e9edf2"
IMBIE_COLOR = "#111111"  # black: no study/category palette uses it
POOLED_COLOR = "#5f6b7a"

TEMPLATE = go.layout.Template(layout=go.Layout(
    font=dict(family=FONT, size=13, color=INK),
    paper_bgcolor="white", plot_bgcolor="white",
    xaxis=dict(gridcolor=GRID, zerolinecolor="#c9d0da", linecolor="#c9d0da", ticks="outside", tickcolor="#c9d0da"),
    yaxis=dict(gridcolor=GRID, zerolinecolor="#c9d0da", linecolor="#c9d0da", ticks="outside", tickcolor="#c9d0da"),
    legend=dict(bgcolor="rgba(255,255,255,0.9)", bordercolor=GRID, borderwidth=1, font=dict(size=12)),
    hoverlabel=dict(font=dict(family=FONT, size=12), bgcolor="white", bordercolor="#c9d0da"),
    margin=dict(l=70, r=20, t=50, b=55),
))

_PUB = RUNS["publication"].to_numpy()
_ICE = RUNS["ice_sheet"].to_numpy()
_HOVER = RUNS["hover"].to_numpy()
_JITTER = np.random.default_rng(42).uniform(-1, 1, len(RUNS))  # fixed per run: points never jump between rebuilds


def rgba(color, alpha):
    c = color.lstrip("#")
    return f"rgba({int(c[0:2], 16)}, {int(c[2:4], 16)}, {int(c[4:6], 16)}, {alpha})"


def _color(dim, cat):
    return POOLED_COLOR if dim is None else COLOR_MAPS[dim].get(cat, POOLED_COLOR)


def units_info(units):
    """(factor, rate label, cumulative label)."""
    if units == "sle":
        return GT_TO_MM_SLE, "mm/yr sea-level equivalent", "mm sea-level equivalent"
    return 1.0, "Gt/yr", "Gt"


def _empty(fig, row, col, text="No simulations selected for this ice sheet"):
    fig.add_annotation(text=text, showarrow=False, font=dict(color=MUTED), row=row, col=col,
                       xref="x domain", yref="y domain", x=0.5, y=0.5)


def _point_style(mask):
    """Big ensembles drawn smaller and fainter so they read as a cloud."""
    sizes, alphas = np.full(len(RUNS), 5.0), np.full(len(RUNS), 0.7)
    for pub in set(_PUB[mask]):
        sel = mask & (_PUB == pub)
        n = sel.sum()
        if n > 30:
            sizes[sel] = 3.5
            alphas[sel] = max(0.15, min(0.7, 20 / n))
    return sizes, alphas


MAX_DRAWN_PER_STUDY = 400


def _drawn(mask):
    """Subset of `mask` actually drawn as dots: at most MAX_DRAWN_PER_STUDY
    evenly spaced runs per study. Statistics (KDEs, medians, boxes, ANOVA)
    always use every run; this only bounds response size for the large
    ensembles (Edwards 2021: 3000 runs per ice sheet)."""
    out = np.zeros(len(RUNS), dtype=bool)
    for pub in set(_PUB[mask]):
        idx = np.flatnonzero(mask & (_PUB == pub))
        if len(idx) > MAX_DRAWN_PER_STUDY:
            idx = idx[np.linspace(0, len(idx) - 1, MAX_DRAWN_PER_STUDY).astype(int)]
        out[idx] = True
    return out


# ── 1. Rates ──────────────────────────────────────────────────────────────

KDE_Y0, KDE_H, STRIP_H, MEDIAN_Y = 0.32, 0.55, 0.16, 0.95


def obs_name(obs):
    return f"{OBS_PRODUCTS[obs]['label']} observed (±2σ)"


def rates_figure(valid, lo, hi, dim, units, show_medians, obs=DEFAULT_OBS, users=()):
    f, rate_label, _ = units_info(units)
    rates = A.run_rates(lo, hi)
    ok = valid & np.isfinite(rates)
    w = A.run_weights(ok)
    sizes, alphas = _point_style(ok)
    drawn = _drawn(ok)
    fig = make_subplots(rows=2, cols=1, vertical_spacing=0.16,
                        subplot_titles=[ICE_SHEET_NAMES[s] for s in ICE_SHEETS])
    legend_seen = set()
    for row, ice in enumerate(ICE_SHEETS, start=1):
        cats = A.category_masks(ok, ice, dim)
        imbie, imbie_se = A.imbie_rate(ice, lo, hi, obs)
        uploads = []  # (upload, finite rates, their simulation names) for this ice sheet
        for u in users:
            if u.ice_sheet == ice:
                r = U.rates(u, lo, hi)
                fin = np.isfinite(r)
                if fin.any():
                    uploads.append((u, r[fin], np.asarray(u.names)[fin]))
        if not cats and not uploads:
            _empty(fig, row, 1)
            continue
        xs_all = np.concatenate([rates[ok & (_ICE == ice)]] + [r for _, r, _ in uploads]) * f
        lo_x, hi_x = min(xs_all.min(), imbie * f), max(xs_all.max(), imbie * f)
        pad = max(0.12 * (hi_x - lo_x), 1e-3)
        grid = np.linspace(lo_x - pad, hi_x + pad, 160)
        for cat, m in cats:
            color = _color(dim, cat)
            x = rates[m] * f
            show = cat not in legend_seen
            legend_seen.add(cat)
            dens = A.kde_curve(x, w[m], grid)
            fig.add_trace(go.Scatter(x=grid, y=np.full_like(grid, KDE_Y0), mode="lines", line=dict(width=0),
                                     hoverinfo="skip", showlegend=False, legendgroup=str(cat)), row=row, col=1)
            fig.add_trace(go.Scatter(
                x=grid, y=KDE_Y0 + dens * KDE_H, mode="lines", fill="tonexty",
                line=dict(color=color, width=1.4), fillcolor=rgba(color, 0.18 if dim else 0.35),
                name=str(cat), legendgroup=str(cat), showlegend=show, hoverinfo="skip"), row=row, col=1)
            d = m & drawn
            fig.add_trace(go.Scatter(
                x=rates[d] * f, y=_JITTER[d] * STRIP_H, mode="markers", legendgroup=str(cat), showlegend=False,
                marker=dict(color=color, size=sizes[d], opacity=alphas[d], line_width=0),
                text=[f"{h}<br><b>Rate: {v:.3g} {rate_label}</b>" for h, v in zip(_HOVER[d], rates[d] * f)],
                hovertemplate="%{text}<extra></extra>"), row=row, col=1)
            if show_medians:
                med = A.weighted_quantiles(x, w[m], [0.5])[0]
                fig.add_trace(go.Scatter(
                    x=[med], y=[MEDIAN_Y], mode="markers", legendgroup=str(cat), showlegend=False,
                    marker=dict(symbol="triangle-down", size=12, color=color, line=dict(width=1, color=INK)),
                    hovertemplate=f"{cat} median: %{{x:.3g}} {rate_label}<extra></extra>"), row=row, col=1)
        for k, (u, r, names) in enumerate(uploads):
            _user_rates(fig, row, u, r * f, names, grid, rate_label, show_medians, k)
        band = sorted([(imbie - 2 * imbie_se) * f, (imbie + 2 * imbie_se) * f])
        fig.add_vrect(x0=band[0], x1=band[1], fillcolor=rgba(IMBIE_COLOR, 0.18), line_width=0, row=row, col=1)
        fig.add_trace(go.Scatter(
            x=[imbie * f] * 2, y=[-0.25, 1.05], mode="lines", line=dict(color=IMBIE_COLOR, width=2.5),
            name=obs_name(obs), legendgroup="imbie", showlegend=row == 1, legendrank=1, hoverinfo="skip"),
            row=row, col=1)
        # The line's own hover points are its ends; this invisible marker puts
        # the IMBIE hover level with the median triangles instead.
        fig.add_trace(go.Scatter(
            x=[imbie * f], y=[MEDIAN_Y], mode="markers", marker=dict(size=14, color="rgba(0,0,0,0)"),
            legendgroup="imbie", showlegend=False,
            hovertemplate=f"{OBS_PRODUCTS[obs]['label']} observed: {imbie * f:.3g} ± {2 * imbie_se * abs(f):.2g} {rate_label}<extra></extra>"),
            row=row, col=1)
        fig.add_vline(x=0, line=dict(color="#9aa3af", dash="dot", width=1), row=row, col=1)
        fig.update_xaxes(title_text=f"Rate of {'sea-level contribution' if units == 'sle' else 'mass change'} ({rate_label})",
                         autorange="reversed" if units == "sle" else True, row=row, col=1)
        fig.update_yaxes(visible=False, range=[-0.3, 1.05], row=row, col=1)
    fig.update_layout(template=TEMPLATE, height=760, uirevision=f"rates-{units}",
                      legend=dict(title=DIM_LABEL.get(dim, "") if dim else None))
    return fig


def _user_rates(fig, row, u, x, names, grid, rate_label, show_medians, k):
    """One of the visitor's uploads in a rates panel: its own density curve,
    dots and median, in its own color, equally weighted (never pooled)."""
    c, label, group = u.color, u.label, f"user:{u.label}"
    dens = A.kde_curve(x, np.ones_like(x), grid)
    fig.add_trace(go.Scatter(x=grid, y=np.full_like(grid, KDE_Y0), mode="lines", line=dict(width=0),
                             hoverinfo="skip", showlegend=False, legendgroup=group), row=row, col=1)
    fig.add_trace(go.Scatter(
        x=grid, y=KDE_Y0 + dens * KDE_H, mode="lines", fill="tonexty", line=dict(color=c, width=2),
        fillcolor=rgba(c, 0.18), name=label, legendgroup=group, showlegend=True, legendrank=2 + k,
        hoverinfo="skip"), row=row, col=1)
    jitter = np.random.default_rng(7 + k).uniform(-1, 1, len(x))
    fig.add_trace(go.Scatter(
        x=x, y=jitter * STRIP_H, mode="markers", legendgroup=group, showlegend=False,
        marker=dict(color=c, size=6 if len(x) <= 30 else 4, opacity=0.8 if len(x) <= 30 else 0.45,
                    line=dict(width=0.5, color="white")),
        text=[f"{label}: {n}<br><b>Rate: {v:.3g} {rate_label}</b>" for n, v in zip(names, x)],
        hovertemplate="%{text}<extra></extra>"), row=row, col=1)
    if show_medians:
        fig.add_trace(go.Scatter(
            x=[np.median(x)], y=[MEDIAN_Y], mode="markers", legendgroup=group, showlegend=False,
            marker=dict(symbol="triangle-down", size=12, color=c, line=dict(width=1, color=INK)),
            hovertemplate=f"{label} median: %{{x:.3g}} {rate_label}<extra></extra>"), row=row, col=1)


# ── 2. Time series ────────────────────────────────────────────────────────

TS_MAX_OUTER_BANDS = 8

def timeseries_figure(valid, lo, hi, dim, units, obs=DEFAULT_OBS, users=()):
    f, _, cum_label = units_info(units)
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.08,
                        subplot_titles=[ICE_SHEET_NAMES[s] for s in ICE_SHEETS])
    legend_seen = set()
    for col, ice in enumerate(ICE_SHEETS, start=1):
        series = A.timeseries(valid, ice, dim, lo, hi)
        for u in users:
            if u.ice_sheet == ice:
                us = U.timeseries(u, lo, hi)
                series = series + ([us] if us else [])
        fig.add_vrect(x0=lo, x1=hi + 1, fillcolor="rgba(250, 204, 21, 0.13)", line_width=0, row=1, col=col)
        ax = "" if col == 1 else str(col)
        fig.add_annotation(text=f"rate window {lo}–{hi}", x=(lo + hi + 1) / 2, y=0.99, xref=f"x{ax}", yref=f"y{ax} domain",
                           yanchor="top", showarrow=False, font=dict(size=11, color="#8a6d00"))
        if not series:
            _empty(fig, 1, col)
        # Layered so lines stay legible with many groups: every group's faint
        # 5-95 % band, then every 25-75 % band, then the medians, then IMBIE
        # on top. Band opacity shrinks as groups are added so overlaps don't
        # pile up into an opaque smear; past TS_MAX_OUTER_BANDS groups the
        # 5-95 % bands are dropped entirely (still given in the hover).
        n = max(len(series), 1)
        a_outer, a_inner = max(0.025, 0.09 / n ** 0.5), max(0.07, 0.24 / n ** 0.5)
        bands = ((1, 3, a_inner),) if n > TS_MAX_OUTER_BANDS else ((0, 4, a_outer), (1, 3, a_inner))
        prepared = []
        for s in series:
            q, ok = s["q"] * f, np.isfinite(s["q"][2])
            prepared.append((s, s.get("color") or _color(dim, s["category"]), s["years"][ok], q[:, ok]))
        for lo_i, hi_i, alpha in bands:
            for s, color, yrs, q in prepared:
                fig.add_trace(go.Scatter(
                    x=np.r_[yrs, yrs[::-1]], y=np.r_[q[hi_i], q[lo_i][::-1]], fill="toself",
                    fillcolor=rgba(color, alpha), line=dict(width=0), hoverinfo="skip",
                    legendgroup=str(s["category"]), showlegend=False), row=1, col=col)
        for s, color, yrs, q in prepared:
            cat = s["category"]
            fig.add_trace(go.Scatter(
                x=yrs, y=q[2], mode="lines", line=dict(color=color, width=2), opacity=0.9,
                name=str(cat), legendgroup=str(cat), showlegend=cat not in legend_seen,
                customdata=np.c_[q[1], q[3], q[0], q[4]],
                hovertemplate=f"<b>{cat}</b> (n={s['n_runs']})<br>%{{x}}: %{{y:.3g}} {cum_label}"
                              f"<br>25–75%: %{{customdata[0]:.3g}} to %{{customdata[1]:.3g}}"
                              f"<br>5–95%: %{{customdata[2]:.3g}} to %{{customdata[3]:.3g}}<extra></extra>"),
                row=1, col=col)
            legend_seen.add(cat)
        ts = A.imbie_timeseries(ice, lo, hi, obs)
        if ts is not None:
            yrs, v, s2 = ts["years"], ts["value"] * f, 2 * ts["sigma"] * abs(f)
            # Hatched (not tinted) uncertainty so it can't be mistaken for a
            # simulation band, edged with thin lines so its extent is crisp.
            fig.add_trace(go.Scatter(
                x=np.r_[yrs, yrs[::-1]], y=np.r_[v + s2, (v - s2)[::-1]], fill="toself",
                fillcolor=rgba(IMBIE_COLOR, 0.05), line=dict(width=0), hoverinfo="skip",
                fillpattern=dict(shape="/", fgcolor=rgba(IMBIE_COLOR, 0.45), size=7, solidity=0.12),
                legendgroup="imbie", showlegend=False), row=1, col=col)
            for edge in (v + s2, v - s2):
                fig.add_trace(go.Scatter(
                    x=yrs, y=edge, mode="lines", line=dict(color=rgba(IMBIE_COLOR, 0.6), width=0.8),
                    hoverinfo="skip", legendgroup="imbie", showlegend=False), row=1, col=col)
            # White halo under the line lifts it off whatever it crosses.
            fig.add_trace(go.Scatter(
                x=yrs, y=v, mode="lines", line=dict(color="white", width=6), hoverinfo="skip",
                legendgroup="imbie", showlegend=False), row=1, col=col)
            fig.add_trace(go.Scatter(
                x=yrs, y=v, mode="lines", line=dict(color=IMBIE_COLOR, width=3),
                name=obs_name(obs), legendgroup="imbie", showlegend=col == 1, legendrank=1,
                hovertemplate=f"{OBS_PRODUCTS[obs]['label']} %{{x:.2f}}: %{{y:.3g}} {cum_label}<extra></extra>"), row=1, col=col)
        fig.add_hline(y=0, line=dict(color="#9aa3af", dash="dot", width=1), row=1, col=col)
        fig.update_xaxes(range=[lo - A.TS_PAD_YEARS, hi + A.TS_PAD_YEARS + 1], title_text="Year", row=1, col=col)
    fig.update_yaxes(title_text=f"Change since {lo} ({cum_label})", row=1, col=1)
    # Vertical legend at the right (as in the rates panel): it scrolls when
    # long, where a horizontal one below the plots squeezes them flat.
    fig.update_layout(template=TEMPLATE, height=460, uirevision=f"ts-{units}",
                      legend=dict(title=DIM_LABEL.get(dim, "") if dim else None))
    return fig


# ── 3. Mass change at 2100 ────────────────────────────────────────────────

# Narrowest a 2100 panel gets: the figure is responsive and its width isn't
# known server-side; at ~1000 px windows (sidebar still beside the content)
# each of the two panels is ~280 px wide. Choosing angles for that width
# keeps labels apart on every wider screen too.
NARROWEST_PANEL_PX = 280


def _tick_style(labels, panel_px=NARROWEST_PANEL_PX):
    """Angle/size for category tick labels so neighbors never overlap: flat
    only if every label fits its slot, else 45 degrees while slanted labels
    stay at least a line apart, else vertical with the font shrunk to the
    slot. Slanted parallel labels are separated by slot * sin(angle),
    regardless of their length."""
    n = max(len(labels), 1)
    slot = panel_px / n
    longest = max((len(str(l)) for l in labels), default=0) * 7.5  # generous px/char at 12 px Inter
    if longest < 0.8 * slot:
        return 0, 12
    if slot * 0.707 >= 18:
        return -45, 12
    return -90, int(max(8, min(12, slot - 2)))


def change_2100_figure(valid, dim, units, users=()):
    f, _, cum_label = units_info(units)
    ok = valid & np.isfinite(A.CHANGE_2100)
    w = A.run_weights(ok)
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.08,
                        subplot_titles=[ICE_SHEET_NAMES[s] for s in ICE_SHEETS])
    label_drop = 0  # px the slanted labels hang below the axis; the figure grows by it
    for col, ice in enumerate(ICE_SHEETS, start=1):
        cats = A.category_masks(ok, ice, dim)
        uploads = []
        for u in users:
            if u.ice_sheet == ice:
                c = U.change_2015_2100(u)
                if np.isfinite(c).any():
                    uploads.append((u, c[np.isfinite(c)]))
        if not cats and not uploads:
            _empty(fig, 1, col)
            continue
        for cat, m in cats:
            x = A.CHANGE_2100[m] * f
            q = A.weighted_quantiles(x, w[m], [0.05, 0.25, 0.5, 0.75, 0.95])
            color = _color(dim, cat)
            fig.add_trace(go.Box(
                x=[str(cat)], q1=[q[1]], median=[q[2]], q3=[q[3]], lowerfence=[q[0]], upperfence=[q[4]],
                fillcolor=rgba(color, 0.35), line=dict(color=color, width=1.6), showlegend=False,
                name=f"{cat} (n={int(m.sum())})", hoverinfo="y+name"), row=1, col=col)
        labels = [c for c, _ in cats]
        for u, uc in uploads:
            q = np.quantile(uc * f, [0.05, 0.25, 0.5, 0.75, 0.95])
            fig.add_trace(go.Box(
                x=[u.label], q1=[q[1]], median=[q[2]], q3=[q[3]], lowerfence=[q[0]], upperfence=[q[4]],
                fillcolor=rgba(u.color, 0.35), line=dict(color=u.color, width=2), showlegend=False,
                name=f"{u.label} (n={len(uc)})", hoverinfo="y+name"), row=1, col=col)
            labels.append(u.label)
        fig.add_hline(y=0, line=dict(color="#9aa3af", dash="dot", width=1), row=1, col=col)
        angle, size = _tick_style(labels)
        longest = max(len(str(c)) for c in labels) * 7 * size / 12
        label_drop = max(label_drop, longest * abs(np.sin(np.radians(angle))))
        fig.update_xaxes(tickangle=angle, tickfont_size=size, automargin=True, row=1, col=col)
    fig.update_yaxes(title_text=f"Change 2015–2100 ({cum_label})", row=1, col=1)
    fig.update_layout(template=TEMPLATE, height=int(400 + label_drop), boxgap=0.35, uirevision=f"y2100-{units}")
    return fig


# ── 4. Bias ───────────────────────────────────────────────────────────────

def bias_figure(b, w, mask, dim, bias_label, obs=DEFAULT_OBS, uploads=()):
    """Weighted box (5/25/50/75/95) + points of bias per category of `dim`,
    sorted by median bias (most negative at top); dashed zero = perfect
    match to IMBIE. `uploads`: [(upload, per-simulation bias)], drawn as
    extra rows at the bottom for comparison only."""
    col = RUNS[dim].to_numpy()
    cats = []
    for c in set(col[mask]):
        m = mask & (col == c)
        cats.append((A.weighted_quantiles(b[m], w[m], [0.5])[0], c, m))
    cats.sort(key=lambda t: t[0])
    sizes, alphas = _point_style(mask)
    drawn = _drawn(mask)
    # Horizontal, one row per category: long labels (study names, the 40-odd
    # GCMs) then sit side by side down the axis instead of colliding along a
    # half-width x-axis, and the figure grows with the number of rows.
    fig = go.Figure()
    for i, (_, c, m) in enumerate(cats):
        color = COLOR_MAPS[dim].get(c, POOLED_COLOR)
        q = A.weighted_quantiles(b[m], w[m], [0.05, 0.25, 0.5, 0.75, 0.95])
        fig.add_trace(go.Box(y=[i], q1=[q[1]], median=[q[2]], q3=[q[3]], lowerfence=[q[0]], upperfence=[q[4]],
                             orientation="h", width=0.6, fillcolor=rgba(color, 0.25), line=dict(color=color, width=1.6),
                             showlegend=False, hoverinfo="x+name", name=f"{c} (n={int(m.sum())})"))
        idx = np.flatnonzero(m & drawn)
        fig.add_trace(go.Scatter(
            x=b[idx], y=i + _JITTER[idx] * 0.26, mode="markers", showlegend=False,
            marker=dict(color=color, size=sizes[idx], opacity=alphas[idx], line_width=0),
            text=[f"{_PUB[j]}: {RUNS.group[j]} / {RUNS.model[j]} / {RUNS.exp[j]}" for j in idx],
            hovertemplate=f"%{{text}}<br>Bias: %{{x:.3g}} {bias_label}<extra></extra>"))
    ticks = [str(c) for _, c, _ in cats]
    for k, (u, ub_all) in enumerate(uploads):
        # Each upload as one extra row at the bottom, for comparison only.
        fin = np.isfinite(ub_all)
        if not fin.any():
            continue
        i, ub, names = len(ticks), ub_all[fin], np.asarray(u.names)[fin]
        q = np.quantile(ub, [0.05, 0.25, 0.5, 0.75, 0.95])
        fig.add_trace(go.Box(y=[i], q1=[q[1]], median=[q[2]], q3=[q[3]], lowerfence=[q[0]], upperfence=[q[4]],
                             orientation="h", width=0.6, fillcolor=rgba(u.color, 0.25),
                             line=dict(color=u.color, width=2), showlegend=False, hoverinfo="x+name",
                             name=f"{u.label} (n={len(ub)})"))
        fig.add_trace(go.Scatter(
            x=ub, y=i + np.random.default_rng(7 + k).uniform(-1, 1, len(ub)) * 0.26, mode="markers", showlegend=False,
            marker=dict(color=u.color, size=5, opacity=0.7, line_width=0),
            text=[f"{u.label}: {n}" for n in names],
            hovertemplate=f"%{{text}}<br>Bias: %{{x:.3g}} {bias_label}<extra></extra>"))
        ticks.append(f"<b>{u.label}</b>")
    fig.add_vline(x=0, line=dict(color=INK, dash="dash", width=1))
    fig.update_layout(template=TEMPLATE, height=max(300, 110 + (34 if len(ticks) <= 12 else 28) * len(ticks)),
                      xaxis_title=f"Bias vs. {OBS_PRODUCTS[obs]['label']} ({bias_label})",
                      yaxis=dict(tickvals=list(range(len(ticks))), ticktext=ticks,
                                 autorange="reversed", automargin=True, showgrid=False, zeroline=False),
                      margin=dict(l=20))
    return fig


def enrichment_figure(rows, dim, frac):
    rows = sorted(rows, key=lambda r: r["share_low_bias"] - r["share_all"])
    cats = [str(r["category"]) for r in rows]
    fig = go.Figure([
        go.Bar(y=cats, x=[100 * r["share_all"] for r in rows], orientation="h", name="All selected runs",
               marker_color="#cdd3dc", hovertemplate="%{y}: %{x:.0f}% of all runs<extra></extra>"),
        go.Bar(y=cats, x=[100 * r["share_low_bias"] for r in rows], orientation="h",
               name=f"Lowest-bias {int(frac * 100)}% of runs", marker_color="#0f9d76",
               hovertemplate="%{y}: %{x:.0f}% of lowest-bias runs<extra></extra>"),
    ])
    fig.update_layout(template=TEMPLATE, barmode="group", height=max(260, 70 + (46 if len(cats) <= 12 else 30) * len(cats)),
                      xaxis_title="Share of runs (weighted, %)", bargap=0.25, bargroupgap=0.05,
                      legend=dict(orientation="h", y=1.08, x=0), yaxis=dict(automargin=True), margin=dict(l=20))
    return fig
