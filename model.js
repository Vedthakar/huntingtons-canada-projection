/* =====================================================================
   Huntington's disease in Canada — age-structured transmission model
   ---------------------------------------------------------------------
   A deterministic, mechanistic cohort model, plus a probabilistic
   sensitivity analysis (Monte Carlo) over its parameter uncertainty.

   Not a statistical fit and not machine learning: every coefficient
   below is a published epidemiological parameter, and the dynamics are
   the Mendelian transmission rule applied year by year to an ageing
   population. The Monte Carlo layer does not learn anything either —
   it propagates stated uncertainty in those parameters through the
   same mechanism, which is how a projection earns an interval.

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

  // Central values for the biological parameters. Each is varied in the
  // Monte Carlo layer over the range given in `priors` below.
  var BIO = {
    piPen: 0.95,        // share of carriers whose allele is fully penetrant (>=40 CAG)
    hdHaz: 0.055,       // excess annual death hazard once symptomatic (~18y median survival)
    tfr: 1.26,          // total fertility rate, Canada 2023
    onsetMedian: 45,
    onsetScale: 8
  };

  /* ---- age schedules --------------------------------------------------
     onset[a]  cumulative probability that a full-penetrance carrier has
               converted to manifest disease by age a (logistic)
     mort[a]   background annual death probability (Gompertz)
     fert[a]   children per person-year, averaged over both sexes
     hazard(a) background mortality plus the HD-specific excess, scaled by
               the share of carriers both fully penetrant and past onset

     Built as a function of the biological parameters so the Monte Carlo
     layer can rebuild them per draw rather than treating them as fixed.
     ------------------------------------------------------------------- */
  function buildSchedules(p) {
    var onset = new Array(A), mort = new Array(A), fert = new Array(A);
    var i, fsum = 0;

    for (i = 0; i < A; i++) {
      onset[i] = i < 15 ? 0 : 1 / (1 + Math.exp(-(i - p.onsetMedian) / p.onsetScale));
      mort[i] = 0.00022 * Math.exp(0.085 * i) + 0.0005;
      if (mort[i] > 0.55) mort[i] = 0.55;
    }
    for (i = 0; i < A; i++) {
      fert[i] = (i >= 15 && i <= 49) ? Math.exp(-Math.pow((i - 30) / 6.2, 2) / 2) : 0;
      fsum += fert[i];
    }
    for (i = 0; i < A; i++) fert[i] = fert[i] / fsum * (p.tfr / 2);

    function hazard(a) {
      return mort[a] + p.piPen * onset[a] * p.hdHaz;
    }

    return { onset: onset, mortality: mort, fertility: fert, hazard: hazard, piPen: p.piPen };
  }

  var SCHED = buildSchedules(BIO);

  /* Stable carrier age distribution under population growth g.
     w[a] is proportional to carrier survivorship to age a, discounted by
     growth so that a growing population has proportionally more young. */
  function stableStructure(g, sched) {
    var s0 = sched || SCHED;
    var w = new Array(A), l = 1, s = 0, a;
    for (a = 0; a < A; a++) {
      w[a] = l * Math.pow(1 + g, -a);
      s += w[a];
      l *= (1 - s0.hazard(a));
      if (l < 0) l = 0;
    }
    for (a = 0; a < A; a++) w[a] /= s;
    return w;
  }

  /* =====================================================================
     run(opts) -> projection
     opts: { cov, eff, start, prev, denovo, growth, bio? }
       cov     0..1    share of at-risk conceptions the policy reaches
       eff     0..1    per-embryo success rate
       start   year     first birth cohort covered
       prev    number   diagnosed cases per 100,000 in 2026
       denovo  0..1    share of incidence arising from new expansions
       growth  0..1    annual population growth rate
       bio     object?  overrides for the biological parameters; when
                        present the age schedules are rebuilt for this run
     ===================================================================== */
  function run(o) {
    var sched = o.bio ? buildSchedules(o.bio) : SCHED;
    var onset = sched.onset, fert = sched.fertility;
    var piPen = sched.piPen;
    var g = o.growth;
    var w = stableStructure(g, sched);
    var a, t;

    // --- calibration ---------------------------------------------------
    // Scale the carrier pool so that the diagnosed prevalence in 2026
    // matches the observed rate.
    var manifestWeight = 0;
    for (a = 0; a < A; a++) manifestWeight += w[a] * piPen * onset[a];
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
        carA += ca[a]; manA += ca[a] * piPen * onset[a];
        carB += cb[a]; manB += cb[a] * piPen * onset[a];
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
        var h = sched.hazard(a - 1);
        ca[a] = ca[a - 1] * (1 - h);
        cb[a] = cb[a - 1] * (1 - h);
      }
      ca[0] = bornA; cb[0] = bornB;
      pop = pop * (1 + g);
    }

    return out;
  }

  /* =====================================================================
     PROBABILISTIC SENSITIVITY ANALYSIS
     ---------------------------------------------------------------------
     A single run of the model is one guess about nine uncertain numbers.
     This layer samples each of them from a stated range, re-runs the whole
     mechanism per draw, and reports percentiles across draws.

     The point is NOT to make the projection more precise. It is to show
     how wide the honest interval is, and which inputs the answer actually
     depends on. A projection quoted without this is a projection quoting
     its own rounding error.

     Ranges below are judgements, stated so they can be argued with:
       prev        BC study point estimate with a ~+/-18% interval
       denovo      literature range for new expansions, 4-10% of incidence
       growth      StatCan low to high immigration scenarios
       eff         editing efficacy: pessimistic to near-perfect
       onsetMedian published medians cluster 40-50
       onsetScale  spread of the onset distribution
       hdHaz       survival post-onset, roughly 15-22 years
       piPen       penetrance of the >=40 CAG threshold
       tfr         fertility, current rate to a modest recovery
     ===================================================================== */
  var priors = {
    prev:        { kind: "normal",   mid: 13.7,  sd: 1.3,   lo: 9,    hi: 19   },
    denovo:      { kind: "uniform",  lo: 0.04,   hi: 0.10 },
    growth:      { kind: "triangle", lo: 0.002,  mid: 0.006, hi: 0.011 },
    eff:         { kind: "triangle", lo: 0.70,   mid: 0.90, hi: 0.99 },
    onsetMedian: { kind: "triangle", lo: 40,     mid: 45,   hi: 50   },
    onsetScale:  { kind: "uniform",  lo: 6.5,    hi: 9.5  },
    hdHaz:       { kind: "uniform",  lo: 0.045,  hi: 0.067 },
    piPen:       { kind: "uniform",  lo: 0.92,   hi: 0.98 },
    tfr:         { kind: "triangle", lo: 1.05,   mid: 1.26, hi: 1.55 }
  };

  // Deterministic RNG (mulberry32) so a reported interval is reproducible.
  function rng(seed) {
    var s = seed >>> 0;
    return function () {
      s = (s + 0x6D2B79F5) >>> 0;
      var x = s;
      x = Math.imul(x ^ (x >>> 15), x | 1);
      x ^= x + Math.imul(x ^ (x >>> 7), x | 61);
      return ((x ^ (x >>> 14)) >>> 0) / 4294967296;
    };
  }

  function sample(p, u1, u2) {
    if (p.kind === "uniform") return p.lo + u1 * (p.hi - p.lo);
    if (p.kind === "triangle") {
      var c = (p.mid - p.lo) / (p.hi - p.lo);
      var v = u1 < c
        ? p.lo + Math.sqrt(u1 * (p.hi - p.lo) * (p.mid - p.lo))
        : p.hi - Math.sqrt((1 - u1) * (p.hi - p.lo) * (p.hi - p.mid));
      return v;
    }
    // normal, via Box-Muller, truncated to [lo, hi]
    var z = Math.sqrt(-2 * Math.log(u1 > 1e-12 ? u1 : 1e-12)) * Math.cos(2 * Math.PI * u2);
    var x = p.mid + z * p.sd;
    if (x < p.lo) x = p.lo;
    if (x > p.hi) x = p.hi;
    return x;
  }

  function pct(sorted, q) {
    if (!sorted.length) return 0;
    var i = (sorted.length - 1) * q;
    var lo = Math.floor(i), hi = Math.ceil(i);
    return lo === hi ? sorted[lo] : sorted[lo] + (sorted[hi] - sorted[lo]) * (i - lo);
  }

  /* montecarlo(opts, n, seed) -> bands + draw-level outputs
       opts  the scenario to test: cov / start are held at the policy
             being argued for; eff and the biology are sampled
       n     number of draws (400 is ample — the bands are stable to
             within a fraction of a percent past ~300)
       seed  integer; the same seed always gives the same interval

     Returns, per year, the 5th / 50th / 95th percentile across draws for
     each of the four curves, plus the distribution of the 2076 reduction
     so the headline number can be quoted as an interval. */
  function montecarlo(o, n, seed) {
    n = n || 400;
    var r = rng(seed == null ? 20260101 : seed);
    var keys = ["aMan", "bMan", "aCar", "bCar"];
    var acc = {}, k, t, d;
    for (k = 0; k < keys.length; k++) {
      acc[keys[k]] = [];
      for (t = 0; t < NY; t++) acc[keys[k]].push(new Array(n));
    }

    var redMan = new Array(n), redCar = new Array(n), averted = new Array(n);
    var sampled = { prev: [], denovo: [], growth: [], eff: [], onsetMedian: [], onsetScale: [], hdHaz: [], piPen: [], tfr: [] };
    var names = Object.keys(sampled);

    for (d = 0; d < n; d++) {
      var draw = {};
      for (k = 0; k < names.length; k++) {
        draw[names[k]] = sample(priors[names[k]], r(), r());
        sampled[names[k]].push(draw[names[k]]);
      }

      var res = run({
        cov: o.cov, start: o.start,
        eff: draw.eff, prev: draw.prev,
        denovo: draw.denovo, growth: draw.growth,
        bio: {
          piPen: draw.piPen, hdHaz: draw.hdHaz, tfr: draw.tfr,
          onsetMedian: draw.onsetMedian, onsetScale: draw.onsetScale
        }
      });

      for (t = 0; t < NY; t++) {
        acc.aMan[t][d] = res.a.man[t];
        acc.bMan[t][d] = res.b.man[t];
        acc.aCar[t][d] = res.a.car[t];
        acc.bCar[t][d] = res.b.car[t];
      }
      var L = NY - 1;
      redMan[d] = res.a.man[L] > 0 ? (res.a.man[L] - res.b.man[L]) / res.a.man[L] : 0;
      redCar[d] = res.a.car[L] > 0 ? (res.a.car[L] - res.b.car[L]) / res.a.car[L] : 0;
      averted[d] = res.a.man[L] - res.b.man[L];
    }

    var bands = {};
    for (k = 0; k < keys.length; k++) {
      var lo = [], mid = [], hi = [];
      for (t = 0; t < NY; t++) {
        var col = acc[keys[k]][t].slice().sort(function (x, y) { return x - y; });
        lo.push(pct(col, 0.05)); mid.push(pct(col, 0.50)); hi.push(pct(col, 0.95));
      }
      bands[keys[k]] = { lo: lo, mid: mid, hi: hi };
    }

    function iv(arr) {
      var s = arr.slice().sort(function (x, y) { return x - y; });
      // min and max are the full observed support, not a percentile. Useful
      // for the strongest honest claim available: "every draw gave at least X".
      return { lo: pct(s, 0.05), mid: pct(s, 0.50), hi: pct(s, 0.95),
               min: s[0], max: s[s.length - 1] };
    }

    /* Which input actually drives the answer?  Rank-correlate each sampled
       parameter against the 2076 carrier reduction. Spearman rather than
       Pearson because several of these act non-linearly. This is a tornado
       plot in numeric form: the top row is the parameter worth arguing about. */
    function ranks(v) {
      var idx = v.map(function (x, i) { return [x, i]; })
                 .sort(function (p1, p2) { return p1[0] - p2[0]; });
      var rk = new Array(v.length);
      for (var i = 0; i < idx.length; i++) rk[idx[i][1]] = i + 1;
      return rk;
    }
    function spearman(x, y) {
      var rx = ranks(x), ry = ranks(y), m = x.length;
      var mx = (m + 1) / 2, num = 0, dx = 0, dy = 0;
      for (var i = 0; i < m; i++) {
        var ax = rx[i] - mx, ay = ry[i] - mx;
        num += ax * ay; dx += ax * ax; dy += ay * ay;
      }
      return (dx > 0 && dy > 0) ? num / Math.sqrt(dx * dy) : 0;
    }

    var drivers = names.map(function (nm) {
      return { param: nm, rho: spearman(sampled[nm], redCar) };
    }).sort(function (p1, p2) { return Math.abs(p2.rho) - Math.abs(p1.rho); });

    var years = [];
    for (t = 0; t < NY; t++) years.push(Y0 + t);

    return {
      n: n, years: years, bands: bands,
      reductionMan2076: iv(redMan),
      reductionCar2076: iv(redCar),
      avertedCases2076: iv(averted),
      drivers: drivers,
      sampled: sampled
    };
  }

  var HDModel = {
    run: run,
    montecarlo: montecarlo,
    priors: priors,
    buildSchedules: buildSchedules,
    schedules: { onset: SCHED.onset, mortality: SCHED.mortality, fertility: SCHED.fertility },
    carrierHazard: SCHED.hazard,
    stableStructure: stableStructure,
    constants: {
      A: A, Y0: Y0, Y1: Y1, NY: NY, POP0: POP0,
      PI_PEN: BIO.piPen, HD_HAZ: BIO.hdHaz, TFR: BIO.tfr,
      ONSET_MEDIAN: BIO.onsetMedian, ONSET_SCALE: BIO.onsetScale
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
