# Ebola (Bundibugyo virus) Outbreak Simulator — Ituri Province, 2026

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://ebola-bdbv-simulator.streamlit.app)

An interactive **Streamlit** decision-support tool that simulates an Ebola disease
outbreak caused by **Bundibugyo virus (BDBV)** in Ituri Province, Democratic Republic
of the Congo. It wraps a compartmental (SEIR-type) transmission model behind public-
health–friendly controls and live visualisations.

> Developed by **The ICMR National Institute of Epidemiology, Chennai** — Indian
> Council of Medical Research, Department of Health Research, Ministry of Health &
> Family Welfare, Government of India.

---

## Features

- **Two-pane layout** — public-health controls on the left, live visualisation on the right.
- **Public-health units** — incubation period, onset→isolation delay, % isolated,
  probability of death (CFR), time-to-death / time-to-recovery, burial delay; the app
  converts these to model rates internally.
- **Honest CFR reporting** — lethality (`p`) is set directly and routed at entry, so the
  realized CFR equals `p` exactly and is **independent of the death/recovery timings**.
  The app plots the real-time **observed** CFR (naive vs resolved) against the true `p`
  to show the right-censoring bias during epidemic growth.
- **Epicurve** with a **Daily / Cumulative** toggle, showing both **cases and deaths**.
- **Two engines** — a fast deterministic RK4 solver (live) and a stochastic
  [`epydemix`](https://pypi.org/project/epydemix/) simulation (averaged over many runs).
- **Tabs** — population groups, observed-vs-true CFR, community-response uptake, and a
  behavioural-response scenario-comparison table.
- **Documentation page** — step-by-step methods, scientific equations alongside plain
  MPH-level English.
- **CSV export** of the full daily trajectory.

## Model (summary)

Eleven compartments: `S, W, A, I, Ii_D, Ii_R, In_D, In_R, D, B, R`. Transmission is
driven by symptomatic (`I`), isolated/non-isolated infectious sub-compartments, and the
deceased-unburied (`D`) — bodies being the dominant BDBV route. Outcome (death vs
recovery) is decided **at entry** via the case-fatality probability `p`; each branch
then drains on its own clock. Transmission rates are auto-calibrated to a target R₀ by
the next-generation method. Full equations, parameter values, ranges and literature
sources are in [`model.md`](model.md).

## Run it in your browser (one click)

Once deployed, click the **Open in Streamlit** badge above to launch the app — no
install required.

**First-time deploy (free, ~1 minute):**

1. Open the prefilled deploy link:
   <https://share.streamlit.io/deploy?repository=drsivambbs/ebola-bdbv-simulator&branch=main&mainModule=app.py>
2. Sign in to Streamlit Community Cloud with the GitHub account that owns this repo.
3. The form is pre-filled (repo `drsivambbs/ebola-bdbv-simulator`, branch `main`,
   file `app.py`). Under **Advanced settings → Custom subdomain**, set it to
   `ebola-bdbv-simulator` so the badge above resolves to
   `https://ebola-bdbv-simulator.streamlit.app`.
4. Click **Deploy**. Streamlit installs `requirements.txt` and serves the app on a
   public URL. Every push to `main` redeploys automatically.

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

Then open the URL Streamlit prints (default http://localhost:8501).
Python 3.10+ recommended. The stochastic engine requires `epydemix`; the default
deterministic engine works without it.

## Disclaimer

This is a research/decision-support tool. Projections are **scenario estimates, not
forecasts**, and must be interpreted alongside field epidemiological assessment.
Several BDBV-specific parameters are borrowed from related ebolaviruses where data are
sparse — see the basis column in `model.md`.
