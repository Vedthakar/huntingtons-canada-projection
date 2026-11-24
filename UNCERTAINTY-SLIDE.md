# The uncertainty slide

Everything here is reproducible with `node verify.js`. All figures: 500 draws,
seed 20260101, broad eligibility (45% coverage, policy from 2035).

---

## The slide

### Title

> **We don't know how many. We do know how much less.**

(Alternative, if the deck's titles are declarative rather than rhetorical:
**"The effect size is better determined than the case count"**)

### The one visual

Two ribbons on a single chart, 2026–2076 — the 90% interval for each scenario.
Live version: the **90% interval** toggle at
https://vedthakar.github.io/huntingtons-canada-projection/

What the eye should catch: the ribbons are **wide** and they **overlap**, yet the
two medians **separate cleanly** after 2035. That picture is the argument. Don't
explain it before they look at it — let the overlap land as a problem, then
resolve it.

### Three numbers, in this order

| | 2076 | |
|---|---|---|
| Baseline cases | **6,266 – 9,771** | the level: ±22% |
| Carrier reduction | **33%** (90% CI 29–38%) | the effect: ±12% |
| Worst of 500 draws | **28%** reduction | the floor |

### The line that does the work

> "Our case count has a ±22% interval, and we'd be lying if we narrowed it.
> Our *effect size* has a ±12% interval — and across 500 runs, not one of them
> came back below a 28% reduction in carriers. So we're not asking you to
> believe a number. We're asking you to notice that the policy conclusion
> doesn't depend on the number."

---

## Why the ratio is tighter than the level (the mechanism — have this ready)

Each of the 500 draws fixes one set of inputs — prevalence, de novo rate,
fertility, onset, everything — and then runs **both scenarios on that same
draw**. So a draw that assumes high prevalence raises the prohibition line *and*
the amendment line together.

When you take the difference as a *proportion*, the shared assumption divides
out. What's left in the interval is only the uncertainty that acts differently
on the two scenarios — essentially just editing efficacy.

This is not a trick of our model. It's the standard reason health-economic
evaluations report incremental effects rather than absolute levels, and it's why
"we don't know the baseline well" is not an objection to a policy comparison.

**Phrase to use:** *shared uncertainty cancels in the ratio.*

---

## The sensitivity ranking — use it defensively

Spearman rank correlation of each sampled input against the 2076 carrier
reduction:

```
efficacy        +0.89
new mutations   -0.33
HD survival     +0.20
pop. growth     +0.18
onset median    -0.16
prevalence      +0.06
fertility       +0.03
penetrance      -0.02
onset spread    -0.00
```

Two things to say about this list:

1. **Only one input matters.** Editing efficacy at ρ = +0.89 is the whole
   argument. That's an honest concession and it's also strategically fine,
   because efficacy is an *empirical* question a clinical framework answers —
   not a values dispute.
2. **Prevalence is at ρ = +0.06, and the true value is zero.** Prevalence scales
   the carrier pool linearly, so it cancels out of the ratio *exactly* — the
   residual 0.06 is Monte Carlo sampling noise, nothing more. So if a judge or
   an opposing team goes after 13.7/100,000, the answer is already on the slide.

> "Double our prevalence estimate and the proportional reduction does not move
> at all — 34.99% either way. The case count doubles; the effect size is
> identical. Prevalence ranks sixth of nine inputs, and its real influence on
> the conclusion is zero."

This is the single most reliable rebuttal in the deck, because it isn't a
judgement call — it's arithmetic. Verify it yourself:

```bash
node -e 'var H=require("./model.js");[13.7,27.4].forEach(function(p){var r=H.run(Object.assign({},H.defaults,H.presets.broad,{prev:p}));console.log(p,((r.a.car[50]-r.b.car[50])/r.a.car[50]*100).toFixed(2)+"%")})'
```

---

## Q&A — prepared answers

**"Where's your confidence interval?"**
On the chart. Nine uncertain inputs, each with a stated prior, resampled 500
times and pushed back through the identical mechanism. 2076 reduction: 33%,
90% CI 29–38%.

**"Isn't this just your assumptions coming back out?"**
Yes — that's exactly what it is, and that's why it's useful. The priors are
written out in `model.js` as nine explicit ranges. You can disagree with a
specific number rather than with "the model." What the analysis adds is showing
you *which* of those nine your disagreement would need to target: realistically,
only efficacy.

**"How do we know 500 runs is enough?"**
Re-run it with a different seed — there's a button. The median moves 0.2%. If it
moved more, we'd need more draws.

**"Is this machine learning / AI?"**
No, and deliberately not. Nothing is learned and no parameter is fitted. This is
*probabilistic sensitivity analysis* — the standard method for propagating
parameter uncertainty in health-policy modelling. Canada has roughly one
reliable Huntington's prevalence estimate and no registry of carrier
conceptions; a model fitted to a single point and extrapolated fifty years
would be the weakest thing in this deck, not the strongest.

**"So the 33% is your real forecast?"**
It's the model's central estimate under our stated priors, at 45% coverage.
Note that 45% coverage is a scenario for answering "where does this lead" — the
amendment we actually propose reaches about 1% of cases and averts roughly 19
diagnosed cases by 2076. The reduction figure is about the mechanism's
behaviour, not a promise about the bill.

---

## Two caveats to state before anyone asks

- **The priors are judgements, not measurements.** The nine ranges were chosen by
  reading the literature, not derived from it. A wider prior would give a wider
  interval. They're published in `model.js` so they're arguable.
- **The draws treat inputs as independent.** Fertility and population growth are
  really linked; we don't model that correlation. It would likely tighten the
  interval slightly, so the version on the slide is the conservative one.

---

## If you only have 30 seconds

> "The case count is uncertain — 6,300 to 9,800 by 2076, and we're not going to
> pretend otherwise. But every one of our 500 runs reduced carriers by at least
> 28%. The level is uncertain; the effect isn't. Judge the policy on the effect."
