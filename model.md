# Disease Transmission Model

A closed-population compartmental model with awareness, isolation, and burial dynamics.

## Compartments

| Symbol | Meaning |
|--------|---------|
| S    | Susceptible |
| W    | Aware — protected; not susceptible and does not transmit |
| A    | Infected — asymptomatic, non-infectious |
| I    | Infected — symptomatic, infectious (initial stage) |
| Iᵢᴰ  | Infected — isolated, **death-bound** |
| Iᵢᴿ  | Infected — isolated, **recovery-bound** |
| Iₙᴰ  | Infected — non-isolated, **death-bound** |
| Iₙᴿ  | Infected — non-isolated, **recovery-bound** |
| D    | Deceased — unburied; transmits until burial |
| B    | Buried — non-infectious, absorbing |
| R    | Recovered — non-infectious |

Transmission is driven by **I, all isolated/non-isolated infectious sub-compartments, and D**.

**Outcome is decided at entry, not by competing rates.** When a case leaves the initial
infectious stage I, it is routed by the *probability of death* **p** (= the CFR, an input) into a
death-bound or a recovery-bound compartment; each then drains on its own clock. This decouples
*how many* die (set by p) from *when* they die/recover (set by the durations), avoiding the
artifact where faster recovery of survivors mechanically lowers the apparent lethality.

## Force of infection

λ = (β₁I + β₂(Iᵢᴰ + Iᵢᴿ) + β₃(Iₙᴰ + Iₙᴿ) + β₄D) · (S / N)

## Awareness (behaviour-driven, endogenous)

ω(t) = min( ω_max, [ κ_c·(I + Iᵢᴰ + Iᵢᴿ + Iₙᴰ + Iₙᴿ) + κ_d·D ] / N )

People adopt protection in proportion to the visible burden — current symptomatic cases and
unburied dead — which for Ebola approximates the recent week's cases and bodies.

## Equations  (γ_d = 1/T_death, γ_r = 1/T_recover)

dS/dt   = − λ − ω(t) S

dW/dt   = ω(t) S

dA/dt   = λ − σ A

dI/dt   = σ A − θ I

dIᵢᴰ/dt = θ · f · p · I − γ_d Iᵢᴰ

dIᵢᴿ/dt = θ · f · (1 − p) · I − γ_r Iᵢᴿ

dIₙᴰ/dt = θ · (1 − f) · p · I − γ_d Iₙᴰ

dIₙᴿ/dt = θ · (1 − f) · (1 − p) · I − γ_r Iₙᴿ

dD/dt   = γ_d (Iᵢᴰ + Iₙᴰ) − b D

dB/dt   = b D

dR/dt   = γ_r (Iᵢᴿ + Iₙᴿ)

Total population is conserved: S + W + A + I + Iᵢᴰ + Iᵢᴿ + Iₙᴰ + Iₙᴿ + D + B + R = N.

By construction the realized CFR equals the input **p**, independent of T_death and T_recover.

## Parameters

Disease: **Ebola disease — Bundibugyo virus (BDBV / *Orthoebolavirus bundibugyoense*)**, the strain of the 2026 DRC–Uganda outbreak (WHO PHEIC declared 17 May 2026). Values are anchored to BDBV-specific literature where it exists (2007–08 Uganda and 2012 DRC outbreaks) and, where BDBV data are absent, to the closest related sources as a last resort — general *Orthoebolavirus* reviews and other Uganda outbreaks (Sudan virus 2022). BDBV-specific gaps are flagged in the basis column. Note two strain-specific facts that shape parameters: **(i) BDBV is markedly less lethal than Zaire (CFR ≈ 33% vs ≈ 67%), and (ii) there is no approved BDBV vaccine or therapeutic**, so awareness (ω), isolation (f), and safe burial (b) are the only control levers. Each rate is a central value with a plausible range and a literature basis; the β's, ω and f are **control/calibration parameters** fixed by fitting to an observed epidemic, with R₀ ≈ 1.5 (range 1.2–2.0) as the primary target.

| Symbol | Meaning | Central | Range | Basis (units: per day unless noted) |
|--------|---------|---------|-------|------|
| N  | Total population | 105,389,731 | — | User-specified (closed population) |
| κ_c | Awareness gain per visible case | 10 | 1–100 | Behavioural response (no lab value); calibrate to observed R(t) decline [6, DON602] |
| κ_d | Awareness gain per unburied death | 50 | 5–500 | Deaths scarier per event; deaths/bodies drive behaviour in BDBV [1] |
| ω_max | Cap on daily protection rate | 0.05 | — | Safety cap; rarely binds at calibrated κ |
| β₁ | Transmission rate, early symptomatic I (community) | 0.15 | 0.08–0.30 | Calibrated so R₀ ≈ 1.5 (BDBV-specific R₀ unpublished; proxy = Uganda outbreaks 1.34–2.7) [2,5] |
| β₂ | Transmission rate, isolated Iᵢ | 0.03 | 0 – 0.07 | ≤ ¼ × β₃ required for effective isolation [8]; target → 0 with barrier nursing |
| β₃ | Transmission rate, non-isolated Iₙ | 0.15 | 0.08–0.30 | ≈ β₁ (same community route); set with β₁ to hit R₀ target |
| β₄ | Transmission rate, deceased D (funeral) | 0.60 | 0.30–1.2 | **Dominant route for BDBV** — handling corpses without protection drove the 2007 outbreak (adj. OR 3.83) [1]; ~2.5 secondary cases per unsafe burial [Tiffany] |
| σ  | Progression A → I (1/incubation) | 0.106 | 0.090–0.167 | No BDBV-specific incubation; proxy = general Ebola 8.5–11.4 d [9,10] and Sudan-Uganda 2022 ≈ 6 d [5] |
| θ  | Exit from initial stage I → Iᵢ/Iₙ (1/onset-to-detection) | 0.22 | 0.14–0.33 | Onset → isolation ≈ 3–5 d [4,5]; note BDBV outbreaks had ~3-month detection lag, so θ is low early [6] |
| f  | Fraction of symptomatic isolated | 0.5 | 0.2–0.8 | Control variable; low at baseline (no vaccine, weak early surveillance), response target > 0.7 [6,8] |
| **p** | **Probability of death (= CFR), set directly** | **0.33** | **0.26–0.40** | BDBV CFR 32.8% (25.8–40.2) [3], 34% / ≈40% in 2007 [1,4]. Input, not an output |
| T_death | Mean time, infectious stage → death (1/γ_d) | 5.0 d | 4–6 | Onset-to-death 9.3 d [9] minus initial stage 1/θ |
| T_recover | Mean time, infectious stage → recovery (1/γ_r) | 8.5 d | 7–10 | Onset-to-recovery 13.0 d [9] minus initial stage 1/θ |
| b  | Burial rate (1/death-to-burial) | 0.4 | 0.33–1.0 | Traditional funeral delay 2–3 d → 0.33–0.5; safe/dignified burial ≤ 1 d → ≥ 1.0 [Tiffany, Checchi] |

### How the parameters are derived from the indicators

Durations come from the consensus delay distributions; the death/recovery split is set by the
**probability p**, not by competing rates:

- σ = 1 / incubation period (A is the latent, non-infectious stage).
- θ = 1 / (symptom-onset → isolation delay) — speed of case-finding; faster detection ⇒ more flow into the isolated branches.
- **p (= CFR) is an input.** A fraction p of cases is routed to the death-bound branch (drains at γ_d = 1/T_death), and (1−p) to the recovery-bound branch (drains at γ_r = 1/T_recover). Realized CFR = p exactly, and is **independent of T_death and T_recover** — so changing care speed changes the *timing* of outcomes, never the lethality.
- b = 1 / (death → burial delay); shortening this is exactly what safe-and-dignified-burial programmes do.

### Reporting CFR honestly (do not present p back as a "result")

Because p is an input, echoing it as a model output is circular. The meaningful, non-circular
CFR outputs are the **real-time observed** ratios, which are biased during a growing epidemic
because deaths lag cases (right-censoring) [Lipsitch 2015; Nishiura 2017; Hauser 2020]:

- **Naive CFR(t) = cumulative deaths / cumulative cases** — *underestimates* p mid-epidemic
  (e.g. reads ≈ 24% when true p = 33%), converging to p only as the wave resolves.
- **Resolved CFR(t) = deaths / (deaths + recoveries)** — much closer to p throughout.

Report both curves against the true p to show the censoring gap; that is the model's genuine
contribution, not the tautological "CFR ≈ 33%".

## Indicators (calibration targets & validation metrics)

These are the empirical quantities the fitted model must reproduce. They are the consensus "indicators" that make the model identifiable rather than free-floating:

| Indicator | Value (BDBV / closest proxy) | Source |
|-----------|-----------------|--------|
| Basic reproduction number R₀ | **No BDBV-specific estimate** (literature gap [9]); proxy = Uganda outbreaks 1.34–2.7, Sudan-Uganda 2022 ≈ 1.25 → target ≈ 1.5 | [2,5] |
| Effective reproduction number R(t) | Tracked over time; must cross 1 as interventions take hold | [7,8] |
| Case-fatality ratio | **32.8% (95% CI 25.8–40.2)** BDBV pooled; 34% / ≈40% / 25% in 2007–08 series | [1,3,4] |
| Incubation period | No BDBV value; general Ebola 8.5–11.4 d, Sudan-Uganda 2022 ≈ 6 d | [5,9,10] |
| Symptom-onset → death | ≈ 9–10 d (general / Sudan-Uganda 2022) | [5,9] |
| Symptom-onset → recovery | ≈ 13 d (general Ebola) | [9] |
| Serial interval | No BDBV value; general Ebola 15.4 d (13.2–17.5) | [9] |
| Dominant transmission route | Handling corpses without protection — adj. OR 3.83 (2007 BDBV) | [1] |
| Secondary cases per unsafe burial | ~2.5 (range 0–20 by district; general Ebola) | [Tiffany] |
| Effect of safe burial | ≥ 40% successful SDB pushes transmission below extinction threshold | [Checchi] |
| Outbreak scale (historical BDBV) | Self-limiting: 2007–08 Uganda 131 cases, 2012 DRC 38 cases | [1,6,8] |

**Calibration workflow:** fix σ, θ, d, ρ, ρ′, b from the duration/CFR indicators above; then fit the β's, ω, and f to the 2026 DRC–Uganda incidence/death time series so the model reproduces R₀ ≈ 1.5 and the early growth (doubling) rate. Use R(t), the BDBV CFR (~33%), and the corpse-dominant transmission pattern as validation. **Caveat:** most durations and R₀ are borrowed from related ebolaviruses because BDBV-specific estimates are sparse — re-fit these as 2026 outbreak data accumulate.

## Sources

BDBV-specific (primary):

1. Wamala JF et al. (2010) *Ebola hemorrhagic fever associated with novel virus strain, Uganda, 2007–2008.* Emerg Infect Dis — BDBV CFR 34%; corpse-handling adj. OR 3.83 (dominant route).
3. Izudi J et al. (2023) *Case fatality rate for Ebola disease, 1976–2022: a meta-analysis.* J Infect Public Health — **BDBV pooled CFR 32.8% (25.8–40.2)**; Zaire 66.6%, Sudan 48.5%.
4. MacNeil A et al. (2010) *Proportion of deaths and clinical features in Bundibugyo Ebola virus infection, Uganda.* Emerg Infect Dis — BDBV CFR ≈ 40%.
6. MacNeil A et al. (2011) *Filovirus outbreak detection and surveillance: lessons from Bundibugyo.* J Infect Dis — 131 cases (2007); ~3-month detection lag.
8b. Hulseberg CE et al. (2021) *Molecular analysis of the 2012 Bundibugyo virus disease outbreak.* Cell Rep Med — 2012 DRC outbreak, 38 confirmed cases.
- Roddy P et al. (2012) *Clinical manifestations and case management of Ebola caused by Bundibugyo, Uganda 2007–08.* PLoS ONE.

Closest-related proxies (used where BDBV data are absent):

2. Muzembo BA et al. (2024) *The basic reproduction number (R₀) of Ebola: systematic review and meta-analysis.* Travel Med Infect Dis — Uganda R₀ range 1.34–2.7.
5. Kabami Z et al. (2024) *Ebola disease outbreak caused by Sudan virus in Uganda, 2022.* Lancet Glob Health — R₀ 1.25, incubation 6 d, onset-to-death 10 d.
7. Thompson RN et al. (2019) *Improved inference of time-varying reproduction numbers...* Epidemics (EpiEstim) — R(t) estimation.
8. Khan A et al. (2015) *Estimating the basic reproductive ratio for the Ebola outbreak in Liberia and Sierra Leone.* Infect Dis Poverty — isolation contact-rate threshold (< ¼).
9. Nash RK et al. (2024) *Ebola virus disease mathematical models and epidemiological parameters: a systematic review.* Lancet Infect Dis — general durations; flags BDBV/non-Zaire data gap.
10. Van Kerkhove MD et al. (2015) *A review of epidemiological parameters from Ebola outbreaks...* Sci Data — incubation, delay distributions.
- Tiffany A et al. (2017) *Estimating the number of secondary Ebola cases from an unsafe burial.* PLoS NTD — ~2.58 secondary cases/unsafe burial.
- Checchi F et al. (2025) *Effect of a safe and dignified burial intervention, DRC 2018–19.* Lancet Glob Health — SDB extinction threshold ~40%.

2026 outbreak context:

- WHO (17 May 2026) *Epidemic of Ebola disease caused by Bundibugyo virus in DRC and Uganda determined a PHEIC.*
- WHO Disease Outbreak News, item **2026-DON602** — Ebola disease caused by Bundibugyo virus, DRC & Uganda.
- CDC HAN / situation summary (June 2026); CDC modelling projecting >20,000 cases in 3 months absent intervention.
