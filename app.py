"""
================================================================================
Streamlit app — Ebola (Bundibugyo virus, BDBV) scenario model, Ituri Province 2026
================================================================================
Interactive wrapper around `ebola_bdbv_ituri_epydemix.py`. The epidemiology is
UNCHANGED — this only exposes the parameters as controls and renders the outputs.

Model logic lives in importable functions (no work at import time):
    calibrate_betas(params)        -> (b1, b2, b3, b4)
    run_deterministic(params)      -> DataFrame (daily trajectory + derived series)
    run_epydemix(params, Nsim)     -> DataFrame (mean trajectory, same columns)
    compute_metrics(df, params)    -> dict (peak_date, peak_new_cases, ...)

Run:
    streamlit run app.py
================================================================================
"""

from __future__ import annotations

import base64
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


# =============================================================================
# DEFAULTS — single source of truth (sourced from model.md / the script)
# =============================================================================
# Controls are stored in PUBLIC-HEALTH units (days, %, counts). The model itself
# runs on per-day rates; to_model_params() converts UI units -> model rates so the
# epidemiology in the run_* functions is completely unchanged.
DEFAULTS = {
    # Population & seeding
    "N": 5_500_000,                 # Ituri Province population (user-specified)
    "A0": 20,                       # initial exposed (incubating) cases
    "I0": 10,                       # initial symptomatic cases
    "start_date": pd.Timestamp("2026-05-17").date(),   # WHO PHEIC declaration
    "horizon_days": 1300,           # projection period (days)
    # Transmission
    "R0_target": 1.5,               # basic reproduction number; proxy = Uganda 1.34-2.7
    "W_COMMUNITY": 1.0,             # relative infectiousness — community case (reference)
    "W_ISOLATED": 0.2,             # relative infectiousness — isolated case
    "W_DECEASED": 4.0,             # relative infectiousness — unsafe burial / body
    # Disease course (public-health units)
    "incubation_days": 9.4,         # 1/sigma  (exposed -> symptomatic)
    "onset_to_isolation_days": 4.5, # 1/theta  (symptom onset -> detection/isolation)
    "pct_isolated": 50.0,           # 100*f    (% of cases isolated)
    # Lethality decoupled from timing: probability of death + two durations.
    "p_death": 0.33,                # p = probability of death (= CFR). Range 0.26-0.40
    "T_death": 5.0,                 # mean time, infectious stage -> death (1/gamma_d)
    "T_recover": 8.5,               # mean time, infectious stage -> recovery (1/gamma_r)
    "death_to_burial_days": 2.5,    # 1/b      (time from death to safe burial)
    # Community / behavioural response
    "kappa_c": 10.0,                # public response to visible cases (strength index)
    "kappa_d": 50.0,                # public response to deaths (strength index)
    "max_adopt_pct": 5.0,           # max % of susceptibles adopting protection per day
    # Engine
    "engine": "Quick estimate (deterministic)",
    "Nsim": 100,
}


# Bump whenever the run_* output schema changes. It rides inside the params dict so
# the st.cache_data key changes too — otherwise edits to helpers like _finish_frame
# (which the cached function doesn't hash) would return stale frames.
_CACHE_VERSION = 2


def to_model_params(ui: dict) -> dict:
    """Convert public-health UI units into the per-day model-rate params dict.

    Durations -> rates (rate = 1/days); percentages -> fractions. Lethality is set
    directly by the probability of death p (= CFR) and routed AT ENTRY, so realized
    CFR = p exactly, independent of the death/recovery timings (T_death, T_recover).
    """
    return {
        "_cache_v": _CACHE_VERSION,
        "N": ui["N"], "A0": ui["A0"], "I0": ui["I0"],
        "start_date": ui["start_date"], "horizon_days": ui["horizon_days"],
        "R0_target": ui["R0_target"],
        "W_COMMUNITY": ui["W_COMMUNITY"], "W_ISOLATED": ui["W_ISOLATED"],
        "W_DECEASED": ui["W_DECEASED"],
        "sigma": 1.0 / ui["incubation_days"],
        "theta": 1.0 / ui["onset_to_isolation_days"],
        "f": ui["pct_isolated"] / 100.0,
        "p": min(max(ui["p_death"], 0.0), 0.99),
        "gamma_d": 1.0 / ui["T_death"],
        "gamma_r": 1.0 / ui["T_recover"],
        "b": 1.0 / ui["death_to_burial_days"],
        "kappa_c": ui["kappa_c"], "kappa_d": ui["kappa_d"],
        "omega_max": ui["max_adopt_pct"] / 100.0,
        "engine": ui["engine"], "Nsim": ui["Nsim"],
    }


# Readable names for the compartment-trajectory legend.
COMPARTMENT_LABELS = {
    "S": "Susceptible", "W": "Protected (behaviour)", "A": "Exposed (incubating)",
    "I": "Symptomatic — early",
    "Ii_D": "Isolated — death-bound", "Ii_R": "Isolated — recovery-bound",
    "In_D": "Non-isolated — death-bound", "In_R": "Non-isolated — recovery-bound",
    "D": "Deceased — unburied", "B": "Buried", "R": "Recovered",
}

COMPARTMENTS = ["S", "W", "A", "I", "Ii_D", "Ii_R", "In_D", "In_R", "D", "B", "R"]


# =============================================================================
# MODEL LOGIC (importable, no side effects)
# =============================================================================
def calibrate_betas(params: dict) -> tuple[float, float, float, float]:
    """Calibrate (b1, b2, b3, b4) to R0_target via the next-generation matrix.

    Next-generation R0 with outcome-at-entry routing. A case spends 1/theta in I
    (beta1), then enters the isolated branch (prob f, beta2) or non-isolated branch
    (prob 1-f, beta3); within each, mean infectious time = p/gamma_d + (1-p)/gamma_r.
    A fraction p eventually reaches D (mean 1/b, beta4):

        R0 = b1/theta + (f*b2 + (1-f)*b3)*(p/gamma_d + (1-p)/gamma_r) + b4*p/b

    Betas keep the same relative structure (b1=b3=W_COMMUNITY*bc, b2=W_ISOLATED*bc,
    b4=W_DECEASED*bc); bc is auto-scaled to hit R0_target.
    """
    wc, wi, wd = params["W_COMMUNITY"], params["W_ISOLATED"], params["W_DECEASED"]
    t_inf = _infectious_time(params)
    coef = (wc * (1.0 / params["theta"])
            + (params["f"] * wi + (1 - params["f"]) * wc) * t_inf
            + wd * params["p"] / params["b"])
    bc = params["R0_target"] / coef          # community transmission base rate
    return (wc * bc,    # b1  (I)
            wi * bc,    # b2  (Ii_D, Ii_R)
            wc * bc,    # b3  (In_D, In_R)
            wd * bc)    # b4  (D)


def _infectious_time(params: dict) -> float:
    """Mean time spent in an isolated/non-isolated branch before D or R."""
    return params["p"] / params["gamma_d"] + (1 - params["p"]) / params["gamma_r"]


def achieved_R0(params: dict, betas: tuple[float, float, float, float]) -> float:
    """Re-derive R0 from the calibrated betas (caption/sanity display)."""
    b1, b2, b3, b4 = betas
    t_inf = _infectious_time(params)
    return (b1 * (1.0 / params["theta"])
            + (params["f"] * b2 + (1 - params["f"]) * b3) * t_inf
            + b4 * params["p"] / params["b"])


def _date_index(params: dict, n_rows: int) -> pd.DatetimeIndex:
    start = pd.Timestamp(params["start_date"])
    return pd.to_datetime([start + pd.Timedelta(days=k) for k in range(n_rows)])


def _finish_frame(Y: np.ndarray, dates, params: dict) -> pd.DataFrame:
    """Assemble the standard daily DataFrame from an (11, T) compartment array."""
    s, w, a, i, iiD, iiR, inD, inR, dd, bb, r = Y
    A0, I0 = params["A0"], params["I0"]
    incidence = params["sigma"] * a                    # new symptomatic cases/day
    new_deaths = params["gamma_d"] * (iiD + inD)       # new deaths/day (inflow to D)
    cum_infected = (s[0] - s - w) + (A0 + I0)          # everyone ever infected
    cum_deaths = dd + bb                               # D + B
    # cumulative symptomatic cases = everyone who has entered I (now downstream of it)
    cum_cases = i + iiD + iiR + inD + inR + dd + bb + r
    # real-time observed CFR (biased low by right-censoring during growth)
    with np.errstate(divide="ignore", invalid="ignore"):
        naive_cfr = np.where(cum_cases > 0, cum_deaths / cum_cases, np.nan)
        resolved_cfr = np.where((cum_deaths + r) > 0,
                                cum_deaths / (cum_deaths + r), np.nan)
    return pd.DataFrame({
        "date": dates,
        "S": s, "W": w, "A": a, "I": i,
        "Ii_D": iiD, "Ii_R": iiR, "In_D": inD, "In_R": inR,
        "D": dd, "B": bb, "R": r,
        "new_cases": incidence,
        "new_deaths": new_deaths,
        "cum_infected": cum_infected,
        "cum_cases": cum_cases,
        "cum_deaths": cum_deaths,
        "cum_recoveries": r,
        "naive_cfr": naive_cfr,
        "resolved_cfr": resolved_cfr,
    })


@st.cache_data(show_spinner=False)
def run_deterministic(params: dict) -> pd.DataFrame:
    """Fast RK4 solver (section 8 of the script), returning a daily DataFrame."""
    b1, b2, b3, b4 = calibrate_betas(params)
    N = params["N"]
    sigma = params["sigma"]; theta = params["theta"]; f = params["f"]
    p = params["p"]; gd = params["gamma_d"]; gr = params["gamma_r"]; b = params["b"]
    kc = params["kappa_c"]; kd = params["kappa_d"]; omega_max = params["omega_max"]
    A0, I0 = params["A0"], params["I0"]
    days = params["horizon_days"]

    def deriv(y):
        s, w, a, i, iiD, iiR, inD, inR, dd, bb, r = y
        omega = min(omega_max, (kc * (i + iiD + iiR + inD + inR) + kd * dd) / N)
        lam = (b1 * i + b2 * (iiD + iiR) + b3 * (inD + inR) + b4 * dd) * s / N
        return np.array([
            -lam - omega * s,                       # dS
            omega * s,                              # dW
            lam - sigma * a,                        # dA
            sigma * a - theta * i,                  # dI
            theta * f * p * i - gd * iiD,           # dIi_D
            theta * f * (1 - p) * i - gr * iiR,     # dIi_R
            theta * (1 - f) * p * i - gd * inD,     # dIn_D
            theta * (1 - f) * (1 - p) * i - gr * inR,  # dIn_R
            gd * (iiD + inD) - b * dd,              # dD
            b * dd,                                 # dB
            gr * (iiR + inR)])                      # dR

    y = np.array([N - A0 - I0, 0, A0, I0, 0, 0, 0, 0, 0, 0, 0], float)
    dt = 0.25
    rows = []
    for step in range(int(days / dt) + 1):
        if abs(step * dt - round(step * dt)) < 1e-9:
            rows.append(y.copy())
        k1 = deriv(y); k2 = deriv(y + dt / 2 * k1)
        k3 = deriv(y + dt / 2 * k2); k4 = deriv(y + dt * k3)
        y = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    Y = np.array(rows).T
    dates = _date_index(params, Y.shape[1])
    return _finish_frame(Y, dates, params)


@st.cache_data(show_spinner=False)
def run_epydemix(params: dict, Nsim: int) -> pd.DataFrame:
    """Stochastic epydemix run, averaged over Nsim simulations.

    Builds the EpiModel exactly as in the script and returns the mean trajectory
    in the same column layout as run_deterministic().
    """
    from epydemix import EpiModel

    b1, b2, b3, b4 = calibrate_betas(params)
    N = params["N"]
    sigma = params["sigma"]; theta = params["theta"]; f = params["f"]
    p = params["p"]; gd = params["gamma_d"]; gr = params["gamma_r"]; b = params["b"]
    kc = params["kappa_c"]; kd = params["kappa_d"]
    A0, I0 = params["A0"], params["I0"]

    parameters = {
        "b1": b1, "b2": b2, "b3": b3, "b4": b4,
        "kappa_c": kc, "kappa_d": kd, "sigma": sigma,
        "theta_f_p": theta * f * p,             "theta_f_1mp": theta * f * (1 - p),
        "theta_1mf_p": theta * (1 - f) * p,     "theta_1mf_1mp": theta * (1 - f) * (1 - p),
        "gamma_d": gd, "gamma_r": gr, "b": b,
    }

    model = EpiModel(
        name="BDBV_Ituri_2026",
        compartments=COMPARTMENTS,
        parameters=parameters,
        use_default_population=True,
        default_population_size=N,
    )

    # Infection: S -> A, mediated by each infectious compartment (forces sum).
    # Both isolated sub-compartments transmit at b2, both non-isolated at b3.
    model.add_transition(source="S", target="A", kind="mediated", params=("b1", "I"))
    model.add_transition(source="S", target="A", kind="mediated", params=("b2", "Ii_D"))
    model.add_transition(source="S", target="A", kind="mediated", params=("b2", "Ii_R"))
    model.add_transition(source="S", target="A", kind="mediated", params=("b3", "In_D"))
    model.add_transition(source="S", target="A", kind="mediated", params=("b3", "In_R"))
    model.add_transition(source="S", target="A", kind="mediated", params=("b4", "D"))
    # Behaviour-driven awareness: S -> W mediated by visible cases & deaths
    model.add_transition(source="S", target="W", kind="mediated", params=("kappa_c", "I"))
    model.add_transition(source="S", target="W", kind="mediated", params=("kappa_c", "Ii_D"))
    model.add_transition(source="S", target="W", kind="mediated", params=("kappa_c", "Ii_R"))
    model.add_transition(source="S", target="W", kind="mediated", params=("kappa_c", "In_D"))
    model.add_transition(source="S", target="W", kind="mediated", params=("kappa_c", "In_R"))
    model.add_transition(source="S", target="W", kind="mediated", params=("kappa_d", "D"))
    # Spontaneous transitions — outcome decided AT ENTRY via p
    model.add_transition(source="A", target="I", kind="spontaneous", params="sigma")
    model.add_transition(source="I", target="Ii_D", kind="spontaneous", params="theta_f_p")
    model.add_transition(source="I", target="Ii_R", kind="spontaneous", params="theta_f_1mp")
    model.add_transition(source="I", target="In_D", kind="spontaneous", params="theta_1mf_p")
    model.add_transition(source="I", target="In_R", kind="spontaneous", params="theta_1mf_1mp")
    model.add_transition(source="Ii_D", target="D", kind="spontaneous", params="gamma_d")
    model.add_transition(source="In_D", target="D", kind="spontaneous", params="gamma_d")
    model.add_transition(source="Ii_R", target="R", kind="spontaneous", params="gamma_r")
    model.add_transition(source="In_R", target="R", kind="spontaneous", params="gamma_r")
    model.add_transition(source="D", target="B", kind="spontaneous", params="b")

    n_groups = len(np.atleast_1d(model.population.Nk))

    def seed(total):
        a = np.zeros(n_groups, dtype=float)
        a[0] = total
        return a

    initial_conditions = {c: seed(0) for c in COMPARTMENTS}
    initial_conditions["A"] = seed(A0)
    initial_conditions["I"] = seed(I0)
    initial_conditions["S"] = seed(N - A0 - I0)

    start_date = pd.Timestamp(params["start_date"]).strftime("%Y-%m-%d")
    end_date = (pd.Timestamp(params["start_date"])
                + pd.Timedelta(days=params["horizon_days"])).strftime("%Y-%m-%d")

    results = model.run_simulations(
        start_date=start_date,
        end_date=end_date,
        initial_conditions_dict=initial_conditions,
        Nsim=Nsim,
        dt=1.0,
    )

    # ---- extract mean compartment trajectories --------------------------------
    # In epydemix 1.2.x each Trajectory.compartments is a dict keyed by
    # "{compartment}_{group}" plus a population-summed "{compartment}_total".
    # We read the "_total" series (summing any age axis) for each compartment.
    dates = pd.to_datetime(np.asarray(results.dates))
    T = len(dates)

    def traj_matrix(tr):
        """Return a (nc, T) array for one trajectory from its compartments dict."""
        comp = tr.compartments
        cols = []
        for c in COMPARTMENTS:
            key = f"{c}_total" if f"{c}_total" in comp else c
            cols.append(np.asarray(comp[key], dtype=float).reshape(-1)[:T])
        return np.vstack(cols)                     # (nc, T)

    stack = np.stack([traj_matrix(tr) for tr in results.trajectories], axis=0)
    Y = stack.mean(axis=0)                         # (nc, T)
    return _finish_frame(Y, dates, params)


def compute_metrics(df: pd.DataFrame, params: dict) -> dict:
    """Headline metrics. CFR is NOT echoed back as a result; instead we report the
    real-time OBSERVED CFR (naive), which under-reads the true p mid-epidemic."""
    incidence = df["new_cases"].to_numpy()
    peak_i = int(np.argmax(incidence))
    peak_date = pd.Timestamp(df["date"].iloc[peak_i])
    peak_new_cases = float(incidence[peak_i])

    total_deaths = float(df["cum_deaths"].iloc[-1])
    cum_infected_end = float(df["cum_infected"].iloc[-1])
    attack_rate = cum_infected_end / params["N"]
    naive_cfr_end = float(df["naive_cfr"].iloc[-1])
    resolved_cfr_end = float(df["resolved_cfr"].iloc[-1])
    # Observed CFR AT THE PEAK — the censored, biased-low value that actually
    # varies with parameters (the end-of-horizon value always converges to p).
    naive_cfr_peak = float(df["naive_cfr"].iloc[peak_i])

    omega_t = (params["kappa_c"] * (df["I"] + df["Ii_D"] + df["Ii_R"]
                                    + df["In_D"] + df["In_R"])
               + params["kappa_d"] * df["D"]) / params["N"]
    omega_peak = float(np.max(omega_t))

    return {
        "peak_date": peak_date,
        "peak_new_cases": peak_new_cases,
        "total_deaths": total_deaths,
        "cum_infected": cum_infected_end,
        "attack_rate": attack_rate,
        "true_cfr": params["p"],
        "naive_cfr_end": naive_cfr_end,
        "naive_cfr_peak": naive_cfr_peak,
        "resolved_cfr_end": resolved_cfr_end,
        "omega_peak": omega_peak,
    }


@st.cache_data(show_spinner=False)
def sensitivity_table(params: dict) -> pd.DataFrame:
    """Behavioural-response sensitivity (section 8): vary kappa_c / kappa_d."""
    grid = [(1, 5), (2, 10), (5, 20), (10, 50), (20, 100), (50, 200)]
    rows = []
    for kc, kd in grid:
        p = dict(params); p["kappa_c"] = float(kc); p["kappa_d"] = float(kd)
        df = run_deterministic(p)
        m = compute_metrics(df, p)
        rows.append({
            "kappa_c": kc, "kappa_d": kd,
            "peak_day": int(np.argmax(df["new_cases"].to_numpy())),
            "peak_new_cases": round(m["peak_new_cases"]),
            "total_deaths": round(m["total_deaths"]),
            "attack_%": round(m["attack_rate"] * 100, 2),
            "infected": round(m["cum_infected"]),
            "obs_CFR_%": round(m["naive_cfr_end"] * 100, 1),
        })
    return pd.DataFrame(rows)


# =============================================================================
# UI HELPERS
# =============================================================================
def _slider(label, key, lo, hi, step, fmt=None, caption=None):
    """Slider that reads its default from DEFAULTS via its key."""
    if key not in st.session_state:
        st.session_state[key] = DEFAULTS[key]
    st.slider(label, lo, hi, step=step, key=key, format=fmt)
    if caption:
        st.caption(caption)


def reset_to_defaults():
    for k in DEFAULTS:
        st.session_state.pop(k, None)


def collect_params() -> dict:
    """Read every control out of session_state into a plain params dict."""
    return {k: st.session_state.get(k, DEFAULTS[k]) for k in DEFAULTS}


# =============================================================================
# BRANDING / THEME
# =============================================================================
# Palette — ICMR / Government-of-India research aesthetic.
NAVY = "#0b2e59"
NAVY_2 = "#13447e"
ACCENT = "#c8102e"        # deep institutional red
INK = "#1b2733"
GRID = "#e6ebf2"

# Compact institutional emblem (generic pathogen motif — not the official
# Government of India State Emblem; used purely as a decorative crest).
_EMBLEM_SVG = """
<svg viewBox="0 0 100 100" width="58" height="58" role="img" aria-label="NIE crest">
  <circle cx="50" cy="50" r="48" fill="#ffffff"/>
  <circle cx="50" cy="50" r="48" fill="none" stroke="#ff9933" stroke-width="3"/>
  <circle cx="50" cy="50" r="44" fill="none" stroke="#138808" stroke-width="3"/>
  <circle cx="50" cy="50" r="38" fill="#0b2e59"/>
  <g stroke="#ffffff" stroke-width="3" stroke-linecap="round">
    <line x1="50" y1="22" x2="50" y2="14"/><line x1="50" y1="78" x2="50" y2="86"/>
    <line x1="22" y1="50" x2="14" y2="50"/><line x1="78" y1="50" x2="86" y2="50"/>
    <line x1="30" y1="30" x2="24" y2="24"/><line x1="70" y1="30" x2="76" y2="24"/>
    <line x1="30" y1="70" x2="24" y2="76"/><line x1="70" y1="70" x2="76" y2="76"/>
  </g>
  <g fill="#ffffff">
    <circle cx="50" cy="14" r="3.4"/><circle cx="50" cy="86" r="3.4"/>
    <circle cx="14" cy="50" r="3.4"/><circle cx="86" cy="50" r="3.4"/>
    <circle cx="24" cy="24" r="3.4"/><circle cx="76" cy="24" r="3.4"/>
    <circle cx="24" cy="76" r="3.4"/><circle cx="76" cy="76" r="3.4"/>
  </g>
  <circle cx="50" cy="50" r="15" fill="#ffffff"/>
  <text x="50" y="55" text-anchor="middle" font-family="Georgia, serif"
        font-size="14" font-weight="700" fill="#0b2e59">NIE</text>
</svg>
"""


@st.cache_data(show_spinner=False)
def _logo_html() -> str:
    """Official ICMR logo as an inline <img> (base64); crest fallback if absent."""
    path = Path(__file__).with_name("icmr_logo.jpeg")
    if path.exists():
        b64 = base64.b64encode(path.read_bytes()).decode("ascii")
        return (f'<div class="icmr-logo">'
                f'<img src="data:image/jpeg;base64,{b64}" alt="ICMR logo"/></div>')
    return f'<div class="icmr-emblem">{_EMBLEM_SVG}</div>'


def inject_css() -> None:
    st.markdown(f"""
    <style>
      .stApp {{ background:#eef2f7; }}
      header[data-testid="stHeader"] {{ background:transparent; height:0; }}
      [data-testid="stToolbar"] {{ right:0.5rem; }}
      .block-container {{ padding-top:1.1rem; padding-bottom:1.5rem; max-width:1500px; }}

      /* ---- masthead ---- */
      .icmr-masthead {{
        display:flex; align-items:center; gap:18px;
        background:linear-gradient(110deg, {NAVY} 0%, {NAVY_2} 100%);
        color:#fff; padding:16px 26px; border-radius:12px 12px 0 0;
        box-shadow:0 4px 16px rgba(11,46,89,.18);
      }}
      .icmr-emblem {{ flex:0 0 auto; line-height:0; }}
      .icmr-logo {{ flex:0 0 auto; background:#fff; padding:7px 13px;
        border-radius:9px; line-height:0; box-shadow:0 1px 4px rgba(0,0,0,.15); }}
      .icmr-logo img {{ height:46px; width:auto; display:block; }}
      .icmr-titles {{ flex:1 1 auto; }}
      .icmr-eyebrow {{ font:600 11px/1 'Segoe UI',system-ui,sans-serif;
        letter-spacing:2.5px; text-transform:uppercase; color:#aac6ec; margin-bottom:5px; }}
      .icmr-org {{ font:700 23px/1.18 Georgia,'Times New Roman',serif; color:#fff; }}
      .icmr-sub {{ font:400 12.5px/1.4 'Segoe UI',system-ui,sans-serif;
        color:#d4e2f4; margin-top:4px; }}
      .icmr-rightmeta {{ flex:0 0 auto; text-align:right; padding-left:14px;
        border-left:1px solid rgba(255,255,255,.18); }}
      .icmr-rt1 {{ font:600 12px/1.3 'Segoe UI',sans-serif; color:#fff; letter-spacing:.4px; }}
      .icmr-rt2 {{ font:400 11px/1.4 'Segoe UI',sans-serif; color:#aac6ec; }}
      .icmr-tricolor {{ height:4px; border-radius:0 0 3px 3px;
        background:linear-gradient(90deg,#ff9933 0 33.3%,#f4f4f4 33.3% 66.6%,#138808 66.6% 100%); }}

      /* ---- app title strip ---- */
      .icmr-apptitle {{ margin:16px 2px 6px; }}
      .icmr-apptitle h1 {{ font:700 27px/1.2 Georgia,serif; color:{NAVY};
        margin:0; padding:0; }}
      .icmr-apptitle p {{ font:400 14px/1.4 'Segoe UI',sans-serif; color:#4a5b6e;
        margin:4px 0 0; }}
      .icmr-pill {{ display:inline-block; font:600 11px/1 'Segoe UI',sans-serif;
        color:{ACCENT}; background:#fdecef; border:1px solid #f5c2cb;
        padding:5px 11px; border-radius:20px; letter-spacing:.4px; margin-top:8px; }}

      /* ---- section headers ---- */
      h2, h3 {{ color:{NAVY} !important; font-family:Georgia,serif !important; }}

      /* ---- metric cards (compact) ---- */
      [data-testid="stMetric"] {{
        background:#fff; border:1px solid #dce4ee; border-left:3px solid {NAVY};
        border-radius:8px; padding:7px 11px;
        box-shadow:0 1px 2px rgba(16,40,73,.05);
      }}
      [data-testid="stMetricLabel"] p {{ font-weight:600; color:#5a6b7e;
        font-size:11px; line-height:1.2; }}
      [data-testid="stMetricValue"] {{ color:{NAVY}; font-weight:700;
        font-size:1.15rem; line-height:1.25; }}

      /* ---- expanders (left controls) ---- */
      [data-testid="stExpander"] {{ border:1px solid #dce4ee; border-radius:10px;
        background:#fff; box-shadow:0 1px 2px rgba(16,40,73,.05); }}
      [data-testid="stExpander"] summary {{ font-weight:600; color:{NAVY}; }}

      /* ---- buttons ---- */
      .stButton button {{ border-radius:8px; font-weight:600; border:1px solid {NAVY};
        color:{NAVY}; background:#fff; }}
      .stButton button:hover {{ background:{NAVY}; color:#fff; border-color:{NAVY}; }}
      [data-testid="stDownloadButton"] button {{ background:{NAVY}; color:#fff;
        border:1px solid {NAVY}; border-radius:8px; font-weight:600; }}
      [data-testid="stDownloadButton"] button:hover {{ background:{NAVY_2}; }}

      /* ---- tabs ---- */
      .stTabs [data-baseweb="tab-list"] {{ gap:4px; }}
      .stTabs [data-baseweb="tab"] {{ font-weight:600; color:#5a6b7e; }}
      .stTabs [aria-selected="true"] {{ color:{NAVY}; }}

      /* ---- documentation page ---- */
      .doc-step-title {{ margin:18px 0 6px; padding:7px 14px; border-radius:8px;
        background:{NAVY}; color:#fff; font:700 15px/1.3 Georgia,serif; }}
      .doc-panehead {{ font:700 13px/1 'Segoe UI',sans-serif; letter-spacing:.5px;
        text-transform:uppercase; color:{NAVY}; padding:4px 2px; }}
      .doc-panehead.sci {{ border-bottom:3px solid {NAVY_2}; }}
      .doc-panehead.plain {{ border-bottom:3px solid {ACCENT}; }}
      .doc-intro {{ font:400 14px/1.6 'Segoe UI',sans-serif; color:#33414f;
        background:#fff; border:1px solid #dce4ee; border-left:4px solid {NAVY};
        border-radius:8px; padding:12px 16px; }}

      /* ---- footer ---- */
      .icmr-footer {{ margin-top:22px; padding:14px 18px; border-top:1px solid #d4ddea;
        font:400 12px/1.6 'Segoe UI',sans-serif; color:#5a6b7e; }}
      .icmr-footer b {{ color:{NAVY}; }}
    </style>
    """, unsafe_allow_html=True)


def render_masthead() -> None:
    st.markdown(f"""
    <div class="icmr-masthead">
      {_logo_html()}
      <div class="icmr-titles">
        <div class="icmr-eyebrow">Developed by</div>
        <div class="icmr-org">The ICMR National Institute of Epidemiology, Chennai</div>
        <div class="icmr-sub">Indian Council of Medical Research &nbsp;·&nbsp;
          Department of Health Research, Ministry of Health &amp; Family Welfare,
          Government of India</div>
      </div>
      <div class="icmr-rightmeta">
        <div class="icmr-rt1">Outbreak Decision-Support</div>
        <div class="icmr-rt2">Compartmental transmission model</div>
      </div>
    </div>
    <div class="icmr-tricolor"></div>
    <div class="icmr-apptitle">
      <h1>Ebola (Bundibugyo virus) Outbreak Simulator</h1>
      <p>Ituri Province, Democratic Republic of the Congo — 2026 scenario projection</p>
      <span class="icmr-pill">BDBV · WHO PHEIC declared 17 May 2026</span>
    </div>
    """, unsafe_allow_html=True)


def style_plot(fig) -> None:
    """Apply the institutional chart theme in place."""
    fig.update_layout(
        font=dict(family="Segoe UI, system-ui, sans-serif", size=13, color=INK),
        title_font=dict(family="Georgia, serif", size=16, color=NAVY),
        paper_bgcolor="#ffffff", plot_bgcolor="#ffffff",
        legend=dict(bgcolor="rgba(255,255,255,0)", borderwidth=0),
        hoverlabel=dict(font_family="Segoe UI, sans-serif"),
        colorway=["#0b2e59", "#13447e", "#2e7da6", "#3aa0a0", "#5fb37a",
                  "#c8102e", "#e07b39", "#8a6db5", "#9aa7b5"],
    )
    fig.update_xaxes(gridcolor=GRID, zeroline=False, linecolor="#c4cedd")
    fig.update_yaxes(gridcolor=GRID, zeroline=False, linecolor="#c4cedd")


# =============================================================================
# DOCUMENTATION PAGE  (scientific | plain-English, step by step)
# =============================================================================
# Each step: (title, [(kind, content), ...] for the science pane, plain-English md).
# kind is "tex" (rendered with st.latex) or "md" (st.markdown).
DOC_STEPS = [
    ("Divide the population into groups", [
        ("md", "The closed population $N$ is partitioned into 11 mutually exclusive "
               "compartments:"),
        ("md", "$S$ susceptible · $W$ protected · $A$ exposed (incubating) · "
               "$I$ symptomatic-early · $I_i^{D},I_i^{R}$ isolated (death- / "
               "recovery-bound) · $I_n^{D},I_n^{R}$ non-isolated · $D$ deceased-"
               "unburied · $B$ buried · $R$ recovered."),
        ("tex", r"S+W+A+I+I_i^{D}+I_i^{R}+I_n^{D}+I_n^{R}+D+B+R = N"),
    ],
     "We track everyone in the population at once, sorting each person into one box "
     "that describes their status today: not yet infected, protecting themselves, "
     "incubating, sick, isolated or not, dead-but-unburied, buried, or recovered. "
     "Nobody leaves the population — they only move between boxes — so the boxes "
     "always add up to the total population."),

    ("Susceptible people catch the virus", [
        ("tex", r"\lambda(t) = \big(\beta_1 I + \beta_2 (I_i^{D}+I_i^{R}) + "
                r"\beta_3 (I_n^{D}+I_n^{R}) + \beta_4 D\big)\,\frac{S}{N}"),
        ("tex", r"\frac{dS}{dt} = -\lambda(t) - \omega(t)\,S "
                r"\qquad \frac{dA}{dt} = \lambda(t) - \sigma A"),
        ("md", "Each infectious group has its own transmissibility $\\beta$. Bodies "
               "($\\beta_4$) are the most infectious route for Bundibugyo virus."),
    ],
     "Healthy (susceptible) people get infected by contact with infectious people. "
     "The chance of infection rises with how many infectious people there are right "
     "now. Not everyone is equally infectious: someone safely isolated spreads very "
     "little, while **unsafe burials and handling of bodies are the biggest driver** "
     "for this virus. Newly infected people move into the 'exposed' box."),

    ("People react and protect themselves", [
        ("tex", r"\omega(t) = \min\!\Big(\omega_{\max},\ "
                r"\frac{\kappa_c\,(I+I_i^{D}+I_i^{R}+I_n^{D}+I_n^{R}) + \kappa_d D}{N}\Big)"),
        ("tex", r"\frac{dW}{dt} = \omega(t)\,S"),
        ("md", "$\\kappa_c,\\kappa_d$ scale how strongly people respond to visible "
               "cases and deaths; $\\omega_{\\max}$ caps the daily uptake."),
    ],
     "As the outbreak becomes visible — more sick people, more funerals — communities "
     "change behaviour: safer burials, avoiding contact, seeking care. We model this "
     "as susceptible people moving into a 'protected' box at a rate that grows with "
     "the burden they can see. Deaths frighten people more than cases, so they drive "
     "a stronger response. There is a ceiling on how fast this can happen."),

    ("Incubation, then symptom onset", [
        ("tex", r"\frac{dA}{dt} = \lambda - \sigma A \qquad "
                r"\frac{dI}{dt} = \sigma A - \theta I"),
        ("md", "$\\sigma = 1/\\text{incubation period}$; "
               "$\\theta = 1/(\\text{onset}\\rightarrow\\text{detection})$."),
    ],
     "After infection there is an incubation period during which the person is **not** "
     "infectious and has no symptoms (the 'exposed' box). On average after about 9 "
     "days symptoms appear and the person becomes infectious. Soon after onset, the "
     "health system either finds and isolates them or it does not — the faster that "
     "happens, the more transmission is cut off."),

    ("Find the case — and decide its outcome up front", [
        ("tex", r"\frac{dI_i^{D}}{dt} = \theta f\,p\,I - \gamma_d I_i^{D} \qquad "
                r"\frac{dI_i^{R}}{dt} = \theta f\,(1-p)\,I - \gamma_r I_i^{R}"),
        ("tex", r"\frac{dI_n^{D}}{dt} = \theta (1\!-\!f)\,p\,I - \gamma_d I_n^{D} \quad "
                r"\frac{dI_n^{R}}{dt} = \theta (1\!-\!f)(1\!-\!p)\,I - \gamma_r I_n^{R}"),
        ("md", "A fraction $f$ are isolated; **independently**, a fraction $p$ are "
               "routed to a death-bound box and $(1-p)$ to a recovery-bound box."),
    ],
     "When a symptomatic case leaves the early stage, two things are decided. "
     "**First**, are they isolated? A share *f* are (the rest stay in the community). "
     "**Second — and this is the key fix — we decide at this moment whether the person "
     "will ultimately die or recover**, using the case-fatality probability *p*. "
     "A share *p* are placed on the 'will die' track and the rest on the 'will "
     "recover' track. Settling the outcome at entry keeps *how many* die separate "
     "from *how long* death or recovery takes."),

    ("Two independent clocks: dying vs recovering", [
        ("tex", r"\gamma_d = \tfrac{1}{T_{\text{death}}}, \quad "
                r"\gamma_r = \tfrac{1}{T_{\text{recover}}}"),
        ("tex", r"\frac{dD}{dt} = \gamma_d (I_i^{D}+I_n^{D}) - bD \qquad "
                r"\frac{dB}{dt} = bD \qquad \frac{dR}{dt} = \gamma_r (I_i^{R}+I_n^{R})"),
    ],
     "People on the 'will die' track die after about *T_death* days and become an "
     "(unburied) body, which keeps transmitting until burial; safe & dignified burial "
     "shortens that window. People on the 'will recover' track recover after about "
     "*T_recover* days. Improving care changes **when** these outcomes happen, never "
     "**how many** people die — that was a flaw in the older design we corrected."),

    ("True vs. observed case-fatality ratio", [
        ("tex", r"\text{CFR}_{\text{true}} = p \quad (\text{by construction})"),
        ("tex", r"\text{CFR}_{\text{naive}}(t) = \frac{D+B}{\text{cumulative cases}} "
                r"\qquad \text{CFR}_{\text{resolved}}(t) = \frac{D+B}{(D+B)+R}"),
    ],
     "Because we set the death probability *p* directly, the **true** case-fatality "
     "ratio is exactly *p*. But during a growing outbreak the **observed** CFR looks "
     "lower than the truth: deaths lag behind newly reported cases (right-censoring), "
     "so dividing deaths by all cases understates lethality until the wave ends. The "
     "‘resolved’ CFR (deaths among people whose outcome is already known) stays much "
     "closer to the truth. The CFR tab shows all three."),

    ("Tune transmission to the target R₀", [
        ("tex", r"R_0 = \frac{\beta_1}{\theta} + \big(f\beta_2 + (1-f)\beta_3\big)"
                r"\left(\frac{p}{\gamma_d} + \frac{1-p}{\gamma_r}\right) + "
                r"\beta_4\,\frac{p}{b}"),
        ("md", "The $\\beta$'s keep fixed relative weights and are scaled together so "
               "this expression equals the chosen $R_0$ (next-generation method)."),
    ],
     "R₀ is the average number of people one case infects in a fully susceptible "
     "population. We don't guess the raw transmission rates; instead we set the R₀ we "
     "want and the model back-calculates the rates that produce it — adding up the "
     "infectiousness contributed during the early stage, during isolation/community "
     "spread, and from bodies. This keeps scenarios comparable when you change other "
     "settings."),
]


def render_documentation() -> None:
    st.markdown("## How the model works")
    st.markdown(
        "<div class='doc-intro'>This is a deterministic compartmental "
        "(SEIR-type) transmission model for Ebola — Bundibugyo virus. Below, each "
        "step appears twice: the <b>scientific formulation</b> on the left and the "
        "<b>same idea in plain language</b> on the right. Full parameter values, "
        "ranges and literature sources are in <i>model.md</i>.</div>",
        unsafe_allow_html=True)
    st.write("")

    hL, hR = st.columns(2, gap="large")
    hL.markdown("<div class='doc-panehead sci'>Scientific formulation</div>",
                unsafe_allow_html=True)
    hR.markdown("<div class='doc-panehead plain'>In plain language · MPH level</div>",
                unsafe_allow_html=True)

    for n, (title, sci_items, plain) in enumerate(DOC_STEPS, start=1):
        st.markdown(f"<div class='doc-step-title'>Step {n} — {title}</div>",
                    unsafe_allow_html=True)
        cL, cR = st.columns(2, gap="large")
        with cL:
            with st.container(border=True):
                for kind, content in sci_items:
                    if kind == "tex":
                        st.latex(content)
                    else:
                        st.markdown(content)
        with cR:
            with st.container(border=True):
                st.markdown(plain)

    st.write("")
    st.markdown(
        "<div class='doc-intro'><b>What the simulator returns:</b> the projected "
        "epicurve (new symptomatic cases/day), peak date and size, total deaths, "
        "attack rate, the population by group over time, community-response uptake, "
        "the observed-vs-true CFR curves, and a scenario-comparison table — all "
        "recomputed live from the parameters on the Simulator page.</div>",
        unsafe_allow_html=True)


# =============================================================================
# PAGE
# =============================================================================
st.set_page_config(
    page_title="BDBV Outbreak Simulator · ICMR-NIE",
    page_icon="🧬",
    layout="wide",
)
inject_css()
render_masthead()

if "nav_page" not in st.session_state:
    st.session_state["nav_page"] = "Simulator"
st.radio("View", ["Simulator", "Documentation"], key="nav_page",
         horizontal=True, label_visibility="collapsed")

if st.session_state["nav_page"] == "Documentation":
    render_documentation()
    st.stop()

left, right = st.columns([1, 2], gap="large")

# ----------------------------------------------------------------------------- LEFT: inputs
with left:
    st.subheader("Parameters")
    st.button("↺ Reset to defaults", on_click=reset_to_defaults,
              width='stretch')

    with st.expander("Population & seeding", expanded=True):
        if "N" not in st.session_state:
            st.session_state["N"] = DEFAULTS["N"]
        st.number_input("Total population at risk", min_value=10_000,
                        max_value=200_000_000, step=100_000, key="N")
        if "A0" not in st.session_state:
            st.session_state["A0"] = DEFAULTS["A0"]
        st.number_input("Initial exposed (incubating) cases", min_value=0,
                        max_value=10_000, step=1, key="A0")
        if "I0" not in st.session_state:
            st.session_state["I0"] = DEFAULTS["I0"]
        st.number_input("Initial symptomatic cases", min_value=0,
                        max_value=10_000, step=1, key="I0")
        if "start_date" not in st.session_state:
            st.session_state["start_date"] = DEFAULTS["start_date"]
        st.date_input("Outbreak start date", key="start_date")
        if "horizon_days" not in st.session_state:
            st.session_state["horizon_days"] = DEFAULTS["horizon_days"]
        st.number_input("Projection period (days)", min_value=100,
                        max_value=3000, step=50, key="horizon_days")

    with st.expander("Transmission", expanded=True):
        _slider("Basic reproduction number (R₀)", "R0_target", 0.8, 3.0, 0.05,
                caption="BDBV R₀ unpublished; proxy = Uganda outbreaks 1.34–2.7 [2,5]")
        _slider("Relative infectiousness — community case", "W_COMMUNITY",
                0.0, 2.0, 0.05, caption="Reference level (= 1.0)")
        _slider("Relative infectiousness — isolated case", "W_ISOLATED",
                0.0, 1.0, 0.05,
                caption="Lower = more effective isolation; ≤ ¼ of community [8]")
        _slider("Relative infectiousness — unsafe burial / body", "W_DECEASED",
                0.0, 8.0, 0.1,
                caption="Dominant BDBV route — unsafe burials (adj. OR 3.83) [1]")

    with st.expander("Disease course (clinical)", expanded=False):
        _slider("Incubation period (days)", "incubation_days", 6.0, 11.5, 0.1,
                "%.1f", "Exposed → symptomatic; general Ebola 8.5–11.4 d [9,10]")
        _slider("Onset → detection / isolation (days)", "onset_to_isolation_days",
                3.0, 7.5, 0.1, "%.1f", "Speed of case-finding [4,5]")
        _slider("Cases isolated (%)", "pct_isolated", 0.0, 100.0, 5.0, "%.0f",
                "Isolation coverage; response target > 70% [6,8]")
        _slider("Probability of death — CFR (p)", "p_death", 0.26, 0.40, 0.01,
                "%.2f", "BDBV ≈ 0.33 (0.26–0.40) [3]. Set directly; realized CFR = p, "
                "independent of the timings below.")
        _slider("Time from onset to death (days)", "T_death", 4.0, 6.0, 0.1, "%.1f",
                "Among fatal cases only; changing it shifts *when*, not *how many* [9]")
        _slider("Time from onset to recovery (days)", "T_recover", 7.0, 10.0, 0.1,
                "%.1f", "Among survivors only; does not affect lethality [9]")
        _slider("Time from death to (safe) burial (days)", "death_to_burial_days",
                1.0, 3.0, 0.1, "%.1f",
                "Safe & dignified burial shortens this [Tiffany, Checchi]")

    with st.expander("Community response", expanded=False):
        _slider("Public response to visible cases", "kappa_c", 0.0, 100.0, 1.0,
                "%.0f", "Higher = faster protective-behaviour uptake per active case")
        _slider("Public response to deaths", "kappa_d", 0.0, 500.0, 5.0, "%.0f",
                "Deaths drive a stronger reaction than cases")
        _slider("Max protective-behaviour adoption (% of susceptibles/day)",
                "max_adopt_pct", 0.0, 20.0, 0.5, "%.1f",
                "Ceiling on how fast people adopt protection (rarely binds)")

    with st.expander("Model mode", expanded=True):
        if "engine" not in st.session_state:
            st.session_state["engine"] = DEFAULTS["engine"]
        st.radio("How to run the projection",
                 ["Quick estimate (deterministic)", "Full simulation (epydemix)"],
                 key="engine",
                 help="Quick estimate recomputes instantly. Full simulation runs "
                      "many stochastic runs and shows an averaged result.")
        if st.session_state["engine"].startswith("Full"):
            if "Nsim" not in st.session_state:
                st.session_state["Nsim"] = DEFAULTS["Nsim"]
            st.slider("Number of simulation runs", 10, 300, step=10, key="Nsim")
            run_stoch = st.button("▶ Run full simulation", width='stretch')
        else:
            run_stoch = False

ui = collect_params()
params = to_model_params(ui)
betas = calibrate_betas(params)

# ----------------------------------------------------------------------------- RIGHT: outputs
with right:
    st.subheader("Projected outbreak — results")
    b1, b2, b3, b4 = betas
    st.caption(
        f"Calibrated transmission rates β₁={b1:.4f}  β₂={b2:.4f}  β₃={b3:.4f}  "
        f"β₄={b4:.4f}  →  achieved R₀ = {achieved_R0(params, betas):.3f}"
    )

    # Always show the fast deterministic result; optionally overlay/replace stochastic.
    det_df = run_deterministic(params)

    stoch_df = None
    is_stochastic = params["engine"].startswith("Full")
    if is_stochastic:
        if run_stoch:
            with st.spinner(f"Running {params['Nsim']} simulations…"):
                try:
                    stoch_df = run_epydemix(params, int(params["Nsim"]))
                    st.session_state["_stoch_cache_key"] = True
                except Exception as e:  # epydemix missing or run error
                    st.error(f"Full simulation failed: {e}")
        else:
            st.info("Full-simulation mode selected — press **Run full simulation** "
                    "in the left panel. Showing the quick estimate until then.")

    # Pick which frame drives the metrics/charts.
    df = stoch_df if stoch_df is not None else det_df
    metrics = compute_metrics(df, params)
    engine_label = ("Full simulation — mean of "
                    f"{int(params['Nsim'])} runs" if stoch_df is not None
                    else "Quick estimate (deterministic)")
    st.caption(f"Engine: **{engine_label}**  ·  Horizon {params['horizon_days']} days "
               f"from {pd.Timestamp(params['start_date']).strftime('%d %b %Y')}")

    # ---- five metric cards ----
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Peak date", metrics["peak_date"].strftime("%Y-%m-%d"))
    c2.metric("New cases on peak day", f"{metrics['peak_new_cases']:,.0f}")
    c3.metric("Total deaths", f"{metrics['total_deaths']:,.0f}")
    c4.metric("Attack rate", f"{metrics['attack_rate'] * 100:.2f}%")
    _gap = (metrics["true_cfr"] - metrics["naive_cfr_peak"]) * 100
    c5.metric("Observed CFR at peak", f"{metrics['naive_cfr_peak'] * 100:.1f}%",
              delta=f"{_gap:+.1f} pts vs true p = {metrics['true_cfr'] * 100:.0f}%",
              delta_color="off",
              help="Naive CFR (cumulative deaths ÷ cumulative cases) on the peak day. "
                   "Right-censoring makes it under-read the true probability of death "
                   "(p) mid-epidemic. It converges to p only after the wave resolves — "
                   "which is why an end-of-horizon value would always just equal p. "
                   "See the CFR tab for the full curves.")

    # ---- epicurve (cases + deaths; daily / cumulative toggle) ----
    view = st.segmented_control(
        "Epicurve view", ["Daily", "Cumulative"], default="Daily",
        key="epicurve_view", label_visibility="collapsed")
    cumulative = (view == "Cumulative")

    fig = go.Figure()
    if cumulative:
        fig.add_trace(go.Scatter(
            x=df["date"], y=df["cum_cases"], mode="lines",
            line=dict(color=ACCENT, width=2.4), fill="tozeroy",
            fillcolor="rgba(200,16,46,.08)", name="Cumulative cases"))
        fig.add_trace(go.Scatter(
            x=df["date"], y=df["cum_deaths"], mode="lines",
            line=dict(color=NAVY, width=2.4), fill="tozeroy",
            fillcolor="rgba(11,46,89,.12)", name="Cumulative deaths"))
        title = "Cumulative cases and deaths"
        ytitle = "Cumulative count"
    else:
        fig.add_trace(go.Scatter(
            x=df["date"], y=df["new_cases"], mode="lines",
            line=dict(color=ACCENT, width=2.4), fill="tozeroy",
            fillcolor="rgba(200,16,46,.08)", name="New cases/day"))
        fig.add_trace(go.Scatter(
            x=df["date"], y=df["new_deaths"], mode="lines",
            line=dict(color=NAVY, width=2.2), name="New deaths/day"))
        fig.add_annotation(
            x=metrics["peak_date"], y=metrics["peak_new_cases"],
            text=f"Peak {metrics['peak_date'].strftime('%d %b %Y')}<br>"
                 f"{metrics['peak_new_cases']:,.0f} cases/day",
            showarrow=True, arrowhead=2, ax=40, ay=-30,
            font=dict(color=NAVY), bgcolor="rgba(255,255,255,.85)")
        title = "Projected epicurve — new cases and deaths per day"
        ytitle = "Per day"
    fig.add_vline(x=metrics["peak_date"], line=dict(color="grey", dash="dash"))
    fig.update_layout(
        title=title, xaxis_title=None, yaxis_title=ytitle,
        height=440, margin=dict(l=10, r=10, t=54, b=64),
        legend=dict(orientation="h", yanchor="top", y=-0.14,
                    x=0.5, xanchor="center"))
    style_plot(fig)
    st.plotly_chart(fig, width='stretch')

    # ---- tabs ----
    tab_comp, tab_cfr, tab_aware, tab_sens = st.tabs(
        ["Population groups", "CFR (observed vs true)",
         "Community response", "Scenario comparison"])

    with tab_comp:
        cfig = go.Figure()
        for c in COMPARTMENTS:
            cfig.add_trace(go.Scatter(x=df["date"], y=df[c], mode="lines",
                                      name=COMPARTMENT_LABELS[c]))
        cfig.update_layout(
            title="Population by group over time", xaxis_title="Date",
            yaxis_title="People", height=460,
            margin=dict(l=10, r=10, t=50, b=10))
        style_plot(cfig)
        st.plotly_chart(cfig, width='stretch')

    with tab_cfr:
        p_true = metrics["true_cfr"]
        ffig = go.Figure()
        ffig.add_trace(go.Scatter(
            x=df["date"], y=[p_true * 100] * len(df), mode="lines",
            line=dict(color=INK, dash="dash", width=1.6),
            name=f"True p = {p_true * 100:.0f}%"))
        ffig.add_trace(go.Scatter(
            x=df["date"], y=df["resolved_cfr"] * 100, mode="lines",
            line=dict(color=NAVY, width=2.4), name="Resolved CFR"))
        ffig.add_trace(go.Scatter(
            x=df["date"], y=df["naive_cfr"] * 100, mode="lines",
            line=dict(color=ACCENT, width=2.4), name="Naive CFR (observed)"))
        ffig.update_layout(
            title="Observed case-fatality ratio over time vs the true value",
            xaxis_title=None, yaxis_title="Case-fatality ratio (%)",
            height=460, margin=dict(l=10, r=10, t=54, b=60),
            legend=dict(orientation="h", yanchor="top", y=-0.14,
                        x=0.5, xanchor="center"))
        style_plot(ffig)
        st.plotly_chart(ffig, width='stretch')
        st.caption("Naive CFR underestimates true lethality mid-epidemic due to "
                   "right-censoring (deaths lag cases); it converges to p only as "
                   "the wave resolves. Resolved CFR stays close to p throughout. "
                   f"At end: naive {metrics['naive_cfr_end']*100:.1f}%, "
                   f"resolved {metrics['resolved_cfr_end']*100:.1f}%, "
                   f"true {p_true*100:.0f}%.")

    with tab_aware:
        omega_t = (params["kappa_c"] * (df["I"] + df["Ii_D"] + df["Ii_R"]
                                        + df["In_D"] + df["In_R"])
                   + params["kappa_d"] * df["D"]) / params["N"]
        omega_t = np.minimum(omega_t, params["omega_max"])
        ofig = go.Figure()
        ofig.add_trace(go.Scatter(x=df["date"], y=omega_t * 100, mode="lines",
                                  line=dict(color=NAVY_2, width=2.4),
                                  fill="tozeroy", fillcolor="rgba(19,68,126,.08)",
                                  name="adoption"))
        ofig.add_hline(y=params["omega_max"] * 100, line=dict(color=ACCENT, dash="dot"),
                       annotation_text="ceiling")
        ofig.update_layout(
            title="Protective-behaviour adoption — rises with visible cases & "
                  "deaths, eases as the wave passes",
            xaxis_title="Date",
            yaxis_title="Susceptibles adopting protection (% per day)",
            height=420, margin=dict(l=10, r=10, t=50, b=10), showlegend=False)
        style_plot(ofig)
        st.plotly_chart(ofig, width='stretch')
        st.caption(f"Peak adoption = {metrics['omega_peak'] * 100:.2f}% of "
                   "susceptibles per day")

    with tab_sens:
        st.caption("How outbreak size changes with the strength of the community "
                   "response, holding all other settings at their current values. "
                   "Stronger response → smaller, later outbreak.")
        sens = sensitivity_table(params).rename(columns={
            "kappa_c": "Response to cases",
            "kappa_d": "Response to deaths",
            "peak_day": "Peak day",
            "peak_new_cases": "Peak new cases",
            "total_deaths": "Total deaths",
            "attack_%": "Attack rate (%)",
            "infected": "Total infected",
            "obs_CFR_%": "Observed CFR (%)",
        })
        st.dataframe(sens, width='stretch', hide_index=True)

    # ---- download ----
    st.download_button(
        "⬇ Download daily trajectory (CSV)",
        data=df.to_csv(index=False).encode("utf-8"),
        file_name="ituri_bdbv_daily.csv",
        mime="text/csv",
        width='stretch')

# ----------------------------------------------------------------------------- footer
st.markdown(f"""
<div class="icmr-footer">
  <b>The ICMR National Institute of Epidemiology, Chennai</b> — Indian Council of
  Medical Research.&nbsp; Bundibugyo virus (BDBV) compartmental transmission model;
  parameters anchored to BDBV-specific literature (see <i>model.md</i>).<br>
  Research decision-support tool — projections are scenario estimates, not forecasts,
  and must be interpreted alongside field epidemiological assessment.
</div>
""", unsafe_allow_html=True)
