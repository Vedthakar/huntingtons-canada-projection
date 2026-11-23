/* =====================================================================
   Huntington's disease in Canada — age-structured transmission model
   ---------------------------------------------------------------------
   A deterministic, mechanistic cohort model. Not a statistical fit and
   not machine learning: every coefficient below is a published
   epidemiological parameter, and the dynamics are the Mendelian
   transmission rule applied year by year to an ageing population.

   Used by index.html (the interactive tool) and how-it-works.html
   (the explainer). Exposes a global `HDModel` and also supports
   CommonJS so the test script can require it.
   ===================================================================== */
(function (root) {
  "use strict";

  // ---- fixed structural constants -----------------------------------
  var A = 101;          // age classes: 0..100
  var Y0 = 2026;
  var Y1 = 2076;
  var NY = Y1 - Y0 + 1;
  var POP0 = 41.5e6;    // Canada, 2026
  var PI_PEN = 0.95;    // share of carriers whose allele is fully penetrant (>=40 CAG)
  var HD_HAZ = 0.055;   // excess annual death hazard once symptomatic (~18y median survival)
  var TFR = 1.26;       // total fertility rate, Canada 2023
  var ONSET_MEDIAN = 45;
  var ONSET_SCALE = 8;

  // ---- age schedules -------------------------------------------------
  // onset[a]  cumulative probability that a full-penetrance carrier has
  //           converted to manifest disease by age a (logistic, median 45)
  // mort[a]   background annual death probability (Gompertz)
  // fert[a]   children per person-year, averaged over both sexes
  var onset = new Array(A);
  var mort = new Array(A);
  var fert = new Array(A);
  var i;

  for (i = 0; i < A; i++) {
    onset[i] = i < 15 ? 0 : 1 / (1 + Math.exp(-(i - ONSET_MEDIAN) / ONSET_SCALE));
    mort[i] = 0.00022 * Math.exp(0.085 * i) + 0.0005;
    if (mort[i] > 0.55) mort[i] = 0.55;
  }

  var fsum = 0;
  for (i = 0; i < A; i++) {
    fert[i] = (i >= 15 && i <= 49) ? Math.exp(-Math.pow((i - 30) / 6.2, 2) / 2) : 0;
    fsum += fert[i];
  }
  for (i = 0; i < A; i++) fert[i] = fert[i] / fsum * (TFR / 2);

  // A carrier dies at the background rate plus an HD-specific excess once
  // symptomatic. The excess is scaled by the share of carriers who are both
  // fully penetrant and past onset at that age.
  function carrierHazard(a) {
    return mort[a] + PI_PEN * onset[a] * HD_HAZ;
  }

  /* Stable carrier age distribution under population growth g.
     w[a] is proportional to carrier survivorship to age a, discounted by
     growth so that a growing population has proportionally more young. */
  function stableStructure(g) {
    var w = new Array(A), l = 1, s = 0, a;
    for (a = 0; a < A; a++) {
      w[a] = l * Math.pow(1 + g, -a);
      s += w[a];
      l *= (1 - carrierHazard(a));
      if (l < 0) l = 0;
    }
    for (a = 0; a < A; a++) w[a] /= s;
    return w;
  }

  /* =====================================================================
     run(opts) -> projection
     opts: { cov, eff, start, prev, denovo, growth }
       cov     0..1   share of at-risk conceptions the policy reaches
       eff     0..1   per-embryo success rate
       start   year    first birth cohort covered
       prev    number  diagnosed cases per 100,000 in 2026
       denovo  0..1   share of incidence arising from new expansions
       growth  0..1   annual population growth rate
     ===================================================================== */
  function run(o) {
    var g = o.growth;
    var w = stableStructure(g);
    var a, t;

    // --- calibration ---------------------------------------------------
    // Scale the carrier pool so that the diagnosed prevalence in 2026
    // matches the observed rate.
    var manifestWeight = 0;
    for (a = 0; a < A; a++) manifestWeight += w[a] * PI_PEN * onset[a];
    var K = (o.prev / 1e5) * POP0 / manifestWeight;   // total carriers, 2026

    // Carrier births needed to hold prevalence-per-capita flat with no
    // intervention. This is what makes the baseline a clean counterfactual:
    // absolute counts rise only because the population rises.
    var required0 = K * w[0] * (1 + g);
    var denovo0 = required0 * o.denovo;
    var inherited0 = required0 - denovo0;

    // Back out the fertility scale that reproduces inherited0 from this
    // age structure, so inherited births respond to the policy thereafter.
    var rawInherited = 0;
    for (a = 0; a < A; a++) rawInherited += K * w[a] * fert[a] * 0.5;
    var fScale = rawInherited > 0 ? inherited0 / rawInherited : 0;
    var denovoPerPop = denovo0 / POP0;   // held constant per capita

    // --- state -----------------------------------------------------------
    var out = {
      years: [], pop: [],
      a: { man: [], car: [] },   // scenario A: prohibition retained
      b: { man: [], car: [] },   // scenario B: germline exception
      averted: [], prevented: [],
      diagnostics: {
        carriers2026: K,
        fertilityScale: fScale,
        denovoBirthsPerYear2026: denovo0,
        inheritedBirthsPerYear2026: inherited0
      }
    };

    var ca = new Array(A), cb = new Array(A);
    for (a = 0; a < A; a++) { ca[a] = K * w[a]; cb[a] = K * w[a]; }
    var pop = POP0, cumPrevented = 0;

    // --- annual loop -------------------------------------------------------
    for (t = 0; t < NY; t++) {
      var yr = Y0 + t;
      var manA = 0, carA = 0, manB = 0, carB = 0;

      for (a = 0; a < A; a++) {
        carA += ca[a]; manA += ca[a] * PI_PEN * onset[a];
        carB += cb[a]; manB += cb[a] * PI_PEN * onset[a];
      }

      out.years.push(yr);
      out.pop.push(pop);
      out.a.man.push(manA); out.a.car.push(carA);
      out.b.man.push(manB); out.b.car.push(carB);
      out.averted.push(manA - manB);
      out.prevented.push(cumPrevented);

      if (t === NY - 1) break;

      // 1. reproduction — each child inherits the expansion with p = 0.5
      var inhA = 0, inhB = 0;
      for (a = 0; a < A; a++) {
        inhA += ca[a] * fert[a] * fScale * 0.5;
        inhB += cb[a] * fert[a] * fScale * 0.5;
      }

      // 2. de novo expansions from intermediate alleles — policy cannot touch these
      var dn = denovoPerPop * pop;

      // 3. the policy filter, applied only to the inherited stream
      var block = (yr >= o.start) ? o.cov * o.eff : 0;
      if (block > 1) block = 1;
      var bornA = inhA + dn;
      var bornB = inhB * (1 - block) + dn;
      cumPrevented += inhB * block;

      // 4. ageing and mortality
      for (a = A - 1; a > 0; a--) {
        var h = carrierHazard(a - 1);
        ca[a] = ca[a - 1] * (1 - h);
        cb[a] = cb[a - 1] * (1 - h);
      }
      ca[0] = bornA; cb[0] = bornB;
      pop = pop * (1 + g);
    }

    return out;
  }

  var HDModel = {
    run: run,
    schedules: { onset: onset, mortality: mort, fertility: fert },
    carrierHazard: carrierHazard,
    stableStructure: stableStructure,
    constants: {
      A: A, Y0: Y0, Y1: Y1, NY: NY, POP0: POP0,
      PI_PEN: PI_PEN, HD_HAZ: HD_HAZ, TFR: TFR,
      ONSET_MEDIAN: ONSET_MEDIAN, ONSET_SCALE: ONSET_SCALE
    },
    defaults: { cov: 0, eff: 0.90, start: 2035, prev: 13.7, denovo: 0.06, growth: 0.006 },
    presets: {
      status:  { cov: 0,    eff: 0.90, start: 2035 },
      narrow:  { cov: 0.01, eff: 0.90, start: 2035 },
      broad:   { cov: 0.45, eff: 0.90, start: 2035 },
      ceiling: { cov: 1,    eff: 1.00, start: 2027 }
    }
  };

  root.HDModel = HDModel;
  if (typeof module !== "undefined" && module.exports) module.exports = HDModel;
})(typeof globalThis !== "undefined" ? globalThis : this);
