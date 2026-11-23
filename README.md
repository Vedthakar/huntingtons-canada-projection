# Huntington's disease in Canada — germline editing policy model

An age-structured transmission model projecting Huntington's disease burden in Canada from **2026 to 2076**, under two scenarios: the current prohibition on germline modification in the *Assisted Human Reproduction Act*, and an amended Act permitting clinical germline gene editing under defined conditions.

Built for a case competition arguing that the AHRA should be amended for autosomal dominant disease where one parent is **homozygous** for the expansion — the case where no unaffected embryo exists and preimplantation genetic testing has nothing to select.

### Open it

| | |
|---|---|
| **The projection tool** | [`index.html`](index.html) — two scenarios, six live parameters, year scrubber |
| **How the model works** | [`how-it-works.html`](how-it-works.html) — diagrams, equations, Q&A prep |

Both run as plain static pages. Clone and open `index.html`, or use the GitHub Pages link in the repo sidebar.

---

## The three findings

**1. The lag is the story.** A conception prevented in 2035 is a case that would have been diagnosed around 2080. Median onset is age 45, so the two diagnosed-case lines stay almost stacked for the first three to four decades and separate only at the far right of the window. Under broad eligibility (45% coverage), 2076 diagnosed cases fall only **11.1%** — but carriers fall **35.0%**. Most of the benefit lands beyond 2076 entirely.

> Lead with this. If a judge finds the flat lines before you explain them it reads as a broken model; if you explain them first it reframes the debate as one about intergenerational duty rather than near-term utility.

**2. The proposed amendment is narrow by design, and the model says so.** Homozygosity is roughly 0.3% of Huntington's cases. At honest coverage (~1%), the amendment averts about **19 diagnosed cases by 2076**. That is not a weakness in the argument — it is the argument. This is a targeted response to absolute medical necessity, not a population health program. The broader presets exist to answer the "where does this lead" question that follows.

**3. It never reaches zero.** New expansions keep arising from intermediate alleles (27–35 CAG) in parents with no family history and no reason to be screened. Even at 100% coverage and 100% efficacy the model floors at **~3,900 diagnosed cases in 2076**, of which 31% of surviving carriers are fresh mutations. Huntington's is reducible, not eliminable. Any claim that germline editing "ends" the disease is false.

---

## This is not a machine learning model

Deliberately. Machine learning finds patterns in training data, and there is no training data here — Canada has roughly one reliable prevalence estimate for Huntington's and no longitudinal registry of carrier conceptions. A model fitted to a single point and extrapolated fifty years would not survive scrutiny.

The generating process is instead **known exactly**: autosomal dominant inheritance is 0.5 per conception, onset follows a measured distribution, populations age at measured rates. When you know the mechanism you encode it rather than infer it. The result is fully auditable — every coefficient traces to a citation, and any assumption can be challenged live by moving a slider.

If asked whether this is AI, the answer is no, and the reason is a strength.

---

## Model structure

Deterministic compartmental cohort model, annual time step, 101 single-year age classes.

Each year, for both scenarios in parallel:

```
1. reproduce      births      = Σ c[a] · f(a) · scale
2. inherit        inherited   = births × 0.5              ← autosomal dominant
3. policy filter  surviving   = inherited × (1 − κ·ε)     ← only if year ≥ start
4. de novo        newborns    = surviving + φ·population  ← policy cannot touch
5. age            c[a] ← c[a−1] · (1 − hazard(a−1));  c[0] ← newborns
```

Diagnosed cases are read off separately through the onset gate: `cases = Σ c[a] · π · F(a)`.

| Symbol | Meaning | Default |
|---|---|---|
| `κ` | share of at-risk conceptions the policy reaches | slider, 0–100% |
| `ε` | per-embryo efficacy | 90% |
| `φ` | de novo share of incidence | 6% |
| `π` | share of carriers with full-penetrance alleles (≥40 CAG) | 0.95 |
| `F(a)` | cumulative onset, logistic | median 45, scale 8 |
| `m(a)` | background mortality, Gompertz | `0.00022·e^0.085a + 0.0005` |
| excess hazard | additional annual death risk once symptomatic | 0.055 (~18y median survival) |
| TFR | total fertility rate | 1.26 |
| growth | annual population growth | 0.60%/yr → 56.0M by 2076 |
| prevalence | diagnosed cases per 100,000 in 2026 | 13.7 |

### Calibration

The carrier pool is scaled so 2026 diagnosed prevalence matches the observed rate, and carrier births are set to the level that holds **prevalence per capita exactly flat** without intervention. This makes the baseline a clean counterfactual: absolute counts rise only because the population rises, so every difference between the two lines is attributable to the policy and nothing else.

Verified: baseline holds 13.70/100k in both 2026 and 2076, with the count rising 5,686 → 7,668 as population grows 41.5M → 56.0M.

---

## Why there is no ancestry term

Prevalence genuinely differs between populations — about 10–14 per 100,000 in European-ancestry populations against under 1 per 100,000 in East Asian populations. That gap is a **founder effect**: differing historical background frequencies of the haplotypes predisposing to CAG expansion.

It is not a difference in transmission. Every carrier passes the expansion to half their conceptions regardless of ancestry. Driving a projection off shifting racial composition would model the wrong mechanism and would be the easiest thing in the deck for an opposing team to attack. The single national prevalence rate, adjustable in the tool, carries the same information without that exposure.

---

## Running the checks

```bash
node verify.js
```

Eleven assertions covering calibration, baseline flatness, the onset lag, the de novo floor, and monotonicity in coverage. All pass on the committed parameters.

## Files

```
index.html         the projection tool
how-it-works.html  explainer with diagrams, equations and Q&A prep
model.js           the model core — shared by both pages and by verify.js
verify.js          numerical sanity checks
```

No build step, no dependencies. `model.js` works in a browser and under Node.

---

## Limits

- No costs, QALYs or budget impact — nothing about cost-effectiveness.
- Coverage is a single number; real uptake would vary by region, income and family awareness, and the hardest-to-reach families are the least likely to have been tested.
- Efficacy is held constant for 50 years. Realistically it improves; flat is the conservative choice.
- Anticipation is ignored — CAG repeats tend to lengthen across generations, especially on paternal transmission, which would shift onset earlier and slightly worsen the no-policy line.
- Reduced-penetrance alleles (36–39 CAG) are folded into a single penetrance constant rather than modelled as their own class.
- The model quantifies one consequence of a policy. It says nothing about whether the policy is right.

State these yourself rather than hoping nobody asks.

---

## Sources

- Prevalence 13.7/100,000 — Fisher & Hayden, British Columbia population study, *Movement Disorders*, 2014.
- Onset distribution, full-penetrance threshold — Caron, Wright & Hayden, *GeneReviews*, 2020.
- Reduced-penetrance alleles 36–39 CAG — Kay et al., *Neurology* 87(3), 2016.
- Homozygosity ~0.3% of cases — Sengupta & Dasgupta, *Clinical Case Reports*, 2026.
- Population base and growth — Statistics Canada population projections, medium-growth scenario.
- Total fertility rate 1.26 — Statistics Canada, 2023.
- De novo expansion from intermediate alleles — standard 5–10% range of incidence.

**Status:** a transparent projection built on published parameters. Not a validated epidemiological forecast, and it should be presented as a model with its assumptions stated.
