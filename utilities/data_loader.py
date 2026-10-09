import netrc
import getpass
import numpy as np
import pandas as pd
import os
import requests
from requests.auth import HTTPDigestAuth
from urllib.parse import urlparse
from urllib.request import HTTPError
from urllib.request import urlopen
from io import BytesIO
from zipfile import ZipFile
from glob import glob
from netCDF4 import Dataset as NC
from pathlib import Path

from .helper import hist_start, hist_end, proj_start, proj_end, proj_time, secpera, ais_exp_dict


def load_imbie_ais():
    """
    Loading the IMBIE Antarctica data set downloaded from
    http://imbie.org/data-files/imbie_dataset-2018_07_23.xlsx
    """
    imbie_df = pd.read_excel(
        "http://imbie.org/data-files/imbie_dataset-2018_07_23.xlsx",
        sheet_name="Antarctica",
        engine="openpyxl",
    )
    imbie = imbie_df[
        [
            "Year",
            "Cumulative ice mass change (Gt)",
            "Cumulative ice mass change uncertainty (Gt)",
        ]
    ].rename(
        columns={
            "Cumulative ice mass change (Gt)": "Cumulative ice sheet mass change (Gt)",
            "Cumulative ice mass change uncertainty (Gt)": "Cumulative ice sheet mass change uncertainty (Gt)",
        }
    )

    for v in [
        "Cumulative ice sheet mass change (Gt)",
    ]:
        imbie[v] -= imbie[imbie["Year"] == proj_start][v].values

    for v in [
        "Cumulative ice sheet mass change uncertainty (Gt)",
    ]:
        imbie[v] -= imbie[v].values[-1]
        imbie[v] *= -1

    imbie["Rate of ice sheet mass change (Gt/yr)"] = np.gradient(
        imbie["Cumulative ice sheet mass change (Gt)"]
    ) / np.gradient(imbie["Year"])
    return imbie


def load_imbie_gis():
    """
    Loading the IMBIE Greenland data set downloaded from
    http://imbie.org/wp-content/uploads/2012/11/imbie_dataset_greenland_dynamics-2020_02_28.xlsx

    """
    imbie_df = pd.read_excel(
        "http://imbie.org/wp-content/uploads/2012/11/imbie_dataset_greenland_dynamics-2020_02_28.xlsx",
        sheet_name="Greenland Ice Mass",
        engine="openpyxl",
    )
    imbie = imbie_df[
        [
            "Year",
            "Cumulative ice sheet mass change (Gt)",
            "Cumulative ice sheet mass change uncertainty (Gt)",
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
    )

    for v in [
        "Cumulative ice sheet mass change (Gt)",
        "Cumulative ice dynamics anomaly (Gt)",
        "Cumulative surface mass balance anomaly (Gt)",
    ]:
        imbie[v] -= imbie[imbie["Year"] == proj_start][v].values

    s = imbie[(imbie["Year"] >= 1980) & (imbie["Year"] < 1990)]
    mass_mean = s["Cumulative ice sheet mass change (Gt)"].mean() / (1990 - 1980)
    smb_mean = s["Cumulative surface mass balance anomaly (Gt)"].mean() / (1990 - 1980)
    imbie[f"Rate of surface mass balance anomaly (Gt/yr)"] += 2 * 1964 / 10
    imbie[f"Rate of ice dynamics anomaly (Gt/yr)"] -= 2 * 1964 / 10
    for v in [
        "Cumulative ice sheet mass change uncertainty (Gt)",
    ]:
        imbie[v] -= imbie[v].values[-1]
        imbie[v] *= -1

    return imbie


def load_ismip6_ais(remove_ctrl=True):
    outpath = "."
    v_dir = "ComputedScalarsPaper"
    url = "https://zenodo.org/record/3940766/files/ComputedScalarsPaper.zip"

    if remove_ctrl:
        ismip6_filename = "ismip6_ais_ctrl_removed.csv.gz"
    else:
        ismip6_filename = "ismip6_ais.csv.gz"
    if os.path.isfile(ismip6_filename):
        df = pd.read_csv(ismip6_filename)
    else:
        print(f"{ismip6_filename} not found locally. Downloading the ISMIP6 archive.")
        # Download only if the archive hasn't already been unpacked here (this
        # used to test for the .zip, which is never saved, so every cache
        # rebuild re-downloaded and re-extracted the whole archive).
        if not os.path.isdir(v_dir) and not os.path.isfile(f"{v_dir}.zip"):
            with urlopen(url) as zipresp:
                with ZipFile(BytesIO(zipresp.read())) as zfile:
                    zfile.extractall(outpath)
        print("   ...and converting to CSV")
        ismip6_ais_to_csv(v_dir, ismip6_filename, remove_ctrl)
        df = pd.read_csv(ismip6_filename)
    return df


def load_ismip6_gis(remove_ctrl=True):
    outpath = "."
    v_dir = "v7_CMIP5_pub"
    url = f"https://zenodo.org/record/3939037/files/{v_dir}.zip"

    if remove_ctrl:
        ismip6_filename = "ismip6_gis_ctrl_removed.csv.gz"
    else:
        ismip6_filename = "ismip6_gis_ctrl.csv.gz"

    if os.path.isfile(ismip6_filename):
        df = pd.read_csv(ismip6_filename)
    else:
        print(f"{ismip6_filename} not found locally. Downloading the ISMIP6 archive.")
        # Download only if the archive hasn't already been unpacked here (this
        # used to test for the .zip, which is never saved, so every cache
        # rebuild re-downloaded and re-extracted the whole archive).
        if not os.path.isdir(v_dir) and not os.path.isfile(f"{v_dir}.zip"):
            with urlopen(url) as zipresp:
                with ZipFile(BytesIO(zipresp.read())) as zfile:
                    zfile.extractall(outpath)
        print("   ...and converting to CSV")
        ismip6_gis_to_csv(v_dir, ismip6_filename, remove_ctrl)
        df = pd.read_csv(ismip6_filename)
    return df


def _gis_mass_af(nc):
    """ISMIP6-Greenland ice mass above flotation (Gt): ivaf (m3) x the
    file's own rhoi (kg/m3) / 1e12."""
    return np.asarray(nc.variables["ivaf"][:], dtype=float) * float(nc.variables["rhoi"][...]) / 1e12


def ismip6_gis_to_csv(basedir, ismip6_filename, remove_ctrl):
    # Now read model output from each of the ISMIP6 files. The information we
    # need is in the file names, not the metadate so this is no fun.
    # Approach is to read each dataset into a dataframe, then concatenate all
    #   dataframes into one Arch dataframe that contains all model runs.
    # Resulting dataframe consists of both historical and projected changes

    ctrl_files = []
    for path in Path(basedir).rglob("*_mm_*_ctrl_proj.nc"):
        ctrl_files.append(path)

    hist_files = []
    for path in Path(basedir).rglob("*_mm_*_historical.nc"):
        hist_files.append(path)

    dfs = []
    for path in Path(basedir).rglob("*_mm_cr_*.nc"):
        # Experiment
        nc = NC(path)
        exp_sle = nc.variables["sle"][:]
        # Mass from ice volume ABOVE FLOTATION x the model's own ice density
        # (rhoi, 900-918 kg/m3 across models), in Gt. Only ice above flotation
        # changes sea level and is what IMBIE's grounded-ice observations
        # track over these timescales; grounded mass (limgr), used here
        # before 2026-10-09, also counts ice below flotation and ran 0.3-22%
        # larger by 2100. ivaf * rhoi equals the files' own limaf (checked to
        # 1e-4), and sle = -limaf / 362.5 Gt per mm.
        exp_mass = _gis_mass_af(nc)
        exp_smb = nc.variables["smb"][:] / 1e12 * secpera

        f = path.name.split(f"scalars_mm_cr_GIS_")[-1].split(".nc")[0].split("_")
        # This is ugly, because of "ITLS_PIK"
        if len(f) == 3:
            group, model, exp = f
        else:
            g1, g2, model, exp = f
            group = f"{g1}_{g2}"

        if exp in ["exp07"]:
            rcp = 26
        else:
            rcp = 85
        # Find the coressponding CTRL Historical simulations
        ctrl_file = [m for m in ctrl_files if (f"_{group}_{model}_" in m.name)][0]
        hist_file = [m for m in hist_files if (f"_{group}_{model}_" in m.name)][0]

        # The last entry of the historical and the first entry of the projection are the same

        # Projection
        nc_ctrl = NC(ctrl_file)
        ctrl_sle = nc_ctrl.variables["sle"][:] - nc_ctrl.variables["sle"][0]
        ctrl_mass = _gis_mass_af(nc_ctrl) - _gis_mass_af(nc_ctrl)[0]
        ctrl_smb = nc_ctrl.variables["smb"][:] / 1e12 * secpera

        # Per email with Heiko on Nov. 13, 2020, stick with just the exp projections alone, without adding back the ctrl projections
        """
        from Heiko:
        "The solution that we chose for ISMIP6 is therefore to remove the ctrl_proj from the projections
        and communicate the numbers as such, i.e. SL contribution for additional forcing after 2014. 
        In our (strong) opinion, the results should never be communicated uncorrected."
        
        Also, point of reference from Goelzer et al., 2020, the ctrl simulations represent mass change
        with the SMB fixed to 1960-1989 levels (no anomaly in SMB) and no change in ice sheet mask.
        So ctrl after the historical spinup represents an abrupt return to an earlier SMB forcing in 2015.
        """

        # Historical
        nc_hist = NC(hist_file)
        hist_sle = nc_hist.variables["sle"][:-1] - nc_hist.variables["sle"][-1]
        hist_mass = _gis_mass_af(nc_hist)[:-1] - _gis_mass_af(nc_hist)[-1]
        hist_smb = nc_hist.variables["smb"][:-1] / 1e12 * secpera
        if remove_ctrl:
            proj_sle = exp_sle
            proj_mass = exp_mass
            proj_smb = exp_smb
        else:
            proj_sle = exp_sle + ctrl_sle
            proj_mass = exp_mass + ctrl_mass
            proj_smb = exp_smb + ctrl_smb

        # Historical simulations start at different years since initialization was left
        # up to the modelers
        hist_time = -np.arange(len(hist_sle))[::-1] + hist_end

        # Let's add the data to the main DataFrame
        m_time = np.hstack((hist_time, proj_time))
        m_sle = -np.hstack((hist_sle, proj_sle)) * 100
        m_sle -= np.interp(proj_start, m_time, m_sle)
        m_mass = np.hstack((hist_mass, proj_mass))
        m_smb = np.cumsum(np.hstack((hist_smb, proj_smb)))
        m_smb -= np.interp(proj_start, m_time, m_smb)
        m_d = m_mass - m_smb
        m_mass_rate = np.gradient(np.hstack((hist_mass, proj_mass)))
        m_smb_rate = np.hstack((hist_smb, proj_smb))
        m_d_rate = m_mass_rate - m_smb_rate
        m_mass -= np.interp(proj_start, m_time, m_mass)

        n = len(m_time)
        dfs.append(
            pd.DataFrame(
                data=np.hstack(
                    [
                        m_time.reshape(-1, 1),
                        m_sle.reshape(-1, 1),
                        m_mass.reshape(-1, 1),
                        m_smb.reshape(-1, 1),
                        m_d.reshape(-1, 1),
                        m_mass_rate.reshape(-1, 1),
                        m_smb_rate.reshape(-1, 1),
                        m_d_rate.reshape(-1, 1),
                        np.repeat(group, n).reshape(-1, 1),
                        np.repeat(model, n).reshape(-1, 1),
                        np.repeat(exp, n).reshape(-1, 1),
                        np.repeat(rcp, n).reshape(-1, 1),
                    ]
                ),
                columns=[
                    "Year",
                    "SLE (cm)",
                    "Cumulative ice sheet mass change (Gt)",
                    "Cumulative surface mass balance anomaly (Gt)",
                    "Cumulative ice dynamics anomaly (Gt)",
                    "Rate of ice sheet mass change (Gt/yr)",
                    "Rate of surface mass balance anomaly (Gt/yr)",
                    "Rate of ice dynamics anomaly (Gt/yr)",
                    "Group",
                    "Model",
                    "Exp",
                    "RCP",
                ],
            )
        )
        # End of working with each model run individually (the path for-loop)

    # Concatenate all DataFrames and convert object types
    df = pd.concat(dfs)
    df = df.astype(
        {
            "Year": float,
            "SLE (cm)": float,
            "Cumulative ice sheet mass change (Gt)": float,
            "Cumulative surface mass balance anomaly (Gt)": float,
            "Cumulative ice dynamics anomaly (Gt)": float,
            "Rate of ice sheet mass change (Gt/yr)": float,
            "Rate of surface mass balance anomaly (Gt/yr)": float,
            "Rate of ice dynamics anomaly (Gt/yr)": float,
            "Model": str,
            "Exp": str,
            "RCP": str,
        }
    )
    df.to_csv(ismip6_filename, compression="gzip")


def _ais_read(path, var):
    """(time, values, rhoi or None) from one ISMIP6-Antarctica computed-scalar file."""
    nc = NC(path)
    rhoi = float(nc.variables["rhoi"][...]) if "rhoi" in nc.variables else None
    return (np.asarray(nc.variables["time"][:], dtype=float), np.asarray(nc.variables[var][:], dtype=float), rhoi)


def ismip6_ais_to_csv(basedir, ismip6_filename, remove_ctrl):
    """One row per (run, year) of ISMIP6-Antarctica (Seroussi et al. 2020)
    historical + projection output, written to `ismip6_filename`.

    Mass is ice volume ABOVE FLOTATION (ivaf) x that model's own ice density
    (rhoi, 900-918 kg/m3 across models; read from each file, or from the same
    model's other files for the few that omit it), in Gt. Before 2026-10-09
    this used total ice volume (ivol) x a fixed 910 kg/m3, which counts
    floating ice shelves -- invisible to sea level and to IMBIE -- and made
    2015-2023 rates ~15x too large (median of exp05 across models: -473 vs
    -32 Gt/yr).

    ivaf is a state at the start of each year (2006.0 ... 2015.0 in hist,
    2016.0 ... 2101.0 in the projections); smb is a flux at mid-year (2015.5,
    ...), so it is interpolated onto the ivaf years rather than stacked as
    separate half-year rows (which also left the mass-rate gradient spanning
    unrelated rows). Rates are d(mass)/dt within each run.

    remove_ctrl: use the *_minus_ctrl_proj projections (ISMIP6's convention)
    with hist rebased to 0 at its end; otherwise the raw projections, with the
    whole series rebased to 0 at proj_start."""
    suffix = "_minus_ctrl_proj" if remove_ctrl else ""
    model_rhoi = {}  # fallback density per (group, model)
    for p in Path(basedir).rglob("computed_ivaf*_AIS_*.nc"):
        r = _ais_read(p, "ivaf")[2]
        if r is not None:
            model_rhoi.setdefault((p.parts[-4], p.parts[-3]), r)

    dfs = []
    for p in sorted(Path(basedir).rglob(f"computed_ivaf{suffix}_AIS_*.nc")):
        if "hist" in p.parent.name or "ctrl" in p.parent.name:
            continue
        group, model, exp = p.parts[-4], p.parts[-3], p.parent.name
        rcp = 26 if exp in ["exp03", "exp07", "expA4", "expA8"] else 85
        t_exp, v_exp, rhoi = _ais_read(p, "ivaf")
        rhoi = rhoi or model_rhoi.get((group, model))
        smb_p = p.with_name(p.name.replace("computed_ivaf", "computed_smb"))
        t_smb, v_smb, _ = _ais_read(smb_p, "smb") if smb_p.exists() else (np.array([]), np.array([]), None)

        hist_dir = Path(basedir) / group / model / f"hist_{ais_exp_dict[exp]}"
        hist_p = hist_dir / f"computed_ivaf_AIS_{group}_{model}_hist_{ais_exp_dict[exp]}.nc"
        if hist_p.exists():
            t_hist, v_hist, _ = _ais_read(hist_p, "ivaf")
            if remove_ctrl:
                v_hist = v_hist - v_hist[-1]
            hist_smb_p = hist_dir / hist_p.name.replace("computed_ivaf", "computed_smb")
            if hist_smb_p.exists():
                th, vh, _ = _ais_read(hist_smb_p, "smb")
                t_smb, v_smb = np.r_[th, t_smb], np.r_[vh, v_smb]
        else:
            t_hist, v_hist = np.array([]), np.array([])

        years = np.r_[t_hist, t_exp]
        mass = np.r_[v_hist, v_exp] * rhoi / 1e12  # m3 above flotation -> Gt
        if not remove_ctrl:
            at0 = np.flatnonzero(years == proj_start)
            mass = mass - (mass[at0[0]] if at0.size else mass[0])
        smb = np.interp(years, t_smb, v_smb) / 1e12 * secpera if t_smb.size else np.full(len(years), np.nan)
        rate = np.gradient(mass, years) if len(years) > 1 else np.full(len(years), np.nan)
        dfs.append(pd.DataFrame({
            "Year": years, "Cumulative ice sheet mass change (Gt)": mass, "Group": group, "Model": model,
            "Exp": exp, "RCP": str(rcp), "Rate of surface mass balance anomaly (Gt/yr)": smb,
            "Rate of ice sheet mass change (Gt/yr)": rate, "Ice density (kg/m3)": rhoi,
        }))
    df = pd.concat(dfs, ignore_index=True)
    df.to_csv(ismip6_filename, compression="gzip")


def plot_historical_partitioning(out_filename, df, imbie):
    def plot_smb(g):
        return ax.plot(
            g[-1]["Year"],
            g[-1]["Rate of surface mass balance anomaly (Gt/yr)"],
            color=simulated_signal_color,
            linewidth=simulated_signal_lw,
            linestyle="solid",
        )

    def plot_d(g):
        return ax.plot(
            g[-1]["Year"],
            g[-1]["Rate of ice dynamics anomaly (Gt/yr)"],
            color=simulated_signal_color,
            linewidth=simulated_signal_lw,
            linestyle="dashed",
        )

    fig = plt.figure()
    ax = fig.add_subplot(111)
    [plot_smb(g) for g in df.groupby(by=["Group", "Model", "Exp"])]
    [plot_d(g) for g in df.groupby(by=["Group", "Model", "Exp"])]

    ax.fill_between(
        imbie["Year"],
        imbie["Rate of surface mass balance anomaly (Gt/yr)"]
        - 1 * imbie["Rate of surface mass balance anomaly uncertainty (Gt/yr)"],
        imbie["Rate of surface mass balance anomaly (Gt/yr)"]
        + 1 * imbie["Rate of surface mass balance anomaly uncertainty (Gt/yr)"],
        color=imbie_sigma_color,
        alpha=0.5,
        linewidth=0,
    )

    ax.fill_between(
        imbie["Year"],
        imbie["Rate of ice dynamics anomaly (Gt/yr)"] - 1 * imbie["Rate of ice dynamics anomaly uncertainty (Gt/yr)"],
        imbie["Rate of ice dynamics anomaly (Gt/yr)"] + 1 * imbie["Rate of ice dynamics anomaly uncertainty (Gt/yr)"],
        color=imbie_sigma_color,
        alpha=0.5,
        linewidth=0,
    )

    ax.plot(
        imbie["Year"],
        imbie["Rate of surface mass balance anomaly (Gt/yr)"],
        color=imbie_signal_color,
        linewidth=imbie_signal_lw,
        linestyle="solid",
    )
    ax.plot(
        imbie["Year"],
        imbie["Rate of ice dynamics anomaly (Gt/yr)"],
        color=imbie_signal_color,
        linewidth=imbie_signal_lw,
        linestyle="dashed",
    )

    ax.axvline(proj_start, color="k", linestyle="dashed", linewidth=grace_signal_lw)
    ax.axhline(0, color="k", linestyle="dotted", linewidth=grace_signal_lw)

    l_smb = mlines.Line2D([], [], color="k", linewidth=0.5, linestyle="solid", label="SMB")
    l_d = mlines.Line2D([], [], color="k", linewidth=0.5, linestyle="dashed", label="D")
    l_imbie = mlines.Line2D(
        [], [], color=imbie_signal_color, linewidth=0.5, linestyle="solid", label="Reconstructed (IMBIE)"
    )
    l_simulated = mlines.Line2D(
        [], [], color=simulated_signal_color, linewidth=0.5, linestyle="solid", label="Simulated"
    )
    legend_1 = ax.legend(
        handles=[l_smb, l_d], loc="upper left", bbox_to_anchor=(0.1, 0.01, 0, 0), bbox_transform=plt.gcf().transFigure
    )
    legend_1.get_frame().set_linewidth(0.0)
    legend_1.get_frame().set_alpha(0.0)

    legend_2 = ax.legend(
        handles=[l_simulated, l_imbie],
        loc="upper left",
        bbox_to_anchor=(0.30, 0.01, 0, 0),
        bbox_transform=plt.gcf().transFigure,
    )
    #    legend_2 = ax.legend(handles=[l_mou19, l_man, l_simulated], loc="upper right")
    legend_2.get_frame().set_linewidth(0.0)
    legend_2.get_frame().set_alpha(0.0)

    # Pylab automacially removes first legend when legend is called a second time.
    # Add legend 1 back
    ax.add_artist(legend_1)

    ax.set_xlabel("Year")
    ax.set_ylabel(f"Flux (Gt/yr)")

    ax.set_xlim(2000, 2014)
    ymin = -750
    ymax = 750
    ax.set_ylim(ymin, ymax)

    set_size(3.2, 2)
    fig.savefig(out_filename, bbox_inches="tight")
