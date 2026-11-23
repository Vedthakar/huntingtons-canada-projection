#!/usr/bin/env node
/* Sanity checks for the projection model.  Run:  node verify.js
   These are the numbers to quote if anyone asks whether the model
   behaves correctly. Every check below must print PASS. */

var HD = require("./model.js");
var base = HD.defaults;
var fails = 0;

function check(name, cond, detail) {
  console.log((cond ? "  PASS  " : "  FAIL  ") + name + (detail ? "   " + detail : ""));
  if (!cond) fails++;
}
function n(x) { return Math.round(x).toLocaleString("en-CA"); }

console.log("\nHuntington's projection model — verification\n" + "-".repeat(58));

// ---- baseline ---------------------------------------------------------
var r = HD.run(base);
var rate0 = r.a.man[0] / r.pop[0] * 1e5;
var rate50 = r.a.man[50] / r.pop[50] * 1e5;

console.log("\nBaseline (prohibition retained)");
console.log("  population      " + (r.pop[0] / 1e6).toFixed(1) + "M -> " + (r.pop[50] / 1e6).toFixed(1) + "M");
console.log("  diagnosed       " + n(r.a.man[0]) + " -> " + n(r.a.man[50]));
console.log("  carriers        " + n(r.a.car[0]) + " -> " + n(r.a.car[50]));
console.log("  rate /100k      " + rate0.toFixed(2) + " -> " + rate50.toFixed(2));
console.log("");

check("calibrates to the observed prevalence", Math.abs(rate0 - base.prev) < 0.01,
      "2026 = " + rate0.toFixed(3) + "/100k");
check("baseline holds per-capita prevalence flat", Math.abs(rate50 - rate0) < 0.02,
      "drift = " + (rate50 - rate0).toFixed(4) + "/100k over 50y");
check("absolute count rises with population", r.a.man[50] > r.a.man[0]);
check("carriers outnumber diagnosed ~4-6x", (r.a.car[0] / r.a.man[0]) > 4 && (r.a.car[0] / r.a.man[0]) < 6,
      "ratio = " + (r.a.car[0] / r.a.man[0]).toFixed(2));

// ---- no policy means no divergence ------------------------------------
check("zero coverage produces identical scenarios",
      Math.abs(r.a.man[50] - r.b.man[50]) < 1e-6);

// ---- the lag ----------------------------------------------------------
var broad = HD.run(Object.assign({}, base, HD.presets.broad));
var gap2045 = (broad.a.man[19] - broad.b.man[19]) / broad.a.man[19];
var gap2076 = (broad.a.man[50] - broad.b.man[50]) / broad.a.man[50];
var carrGap2076 = (broad.a.car[50] - broad.b.car[50]) / broad.a.car[50];

console.log("\nBroad eligibility (45% coverage, policy from 2035)");
console.log("  2076 diagnosed  " + n(broad.a.man[50]) + " vs " + n(broad.b.man[50]) +
            "   (" + (gap2076 * 100).toFixed(1) + "% lower)");
console.log("  2076 carriers   " + n(broad.a.car[50]) + " vs " + n(broad.b.car[50]) +
            "   (" + (carrGap2076 * 100).toFixed(1) + "% lower)");
console.log("  births averted  " + n(broad.prevented[50]) + " cumulative");
console.log("");

check("onset lag suppresses the near-term effect", gap2045 < 0.02,
      "2045 gap = " + (gap2045 * 100).toFixed(2) + "%");
check("carrier effect far exceeds diagnosed effect by 2076", carrGap2076 > gap2076 * 2,
      "carriers " + (carrGap2076 * 100).toFixed(1) + "% vs cases " + (gap2076 * 100).toFixed(1) + "%");

// ---- the de novo floor ------------------------------------------------
var ceil = HD.run(Object.assign({}, base, HD.presets.ceiling));
console.log("\nTheoretical ceiling (100% coverage, 100% efficacy, from 2027)");
console.log("  2076 diagnosed  " + n(ceil.a.man[50]) + " vs " + n(ceil.b.man[50]));
console.log("  2076 carriers   " + n(ceil.a.car[50]) + " vs " + n(ceil.b.car[50]));
console.log("");

check("even the ceiling cannot reach zero", ceil.b.man[50] > 100,
      "floor = " + n(ceil.b.man[50]) + " cases held up by new mutations");

// Under the ceiling the carriers left in 2076 are of two kinds: people
// already alive in 2026 who simply have not died yet, and fresh expansions.
// Turning de novo off isolates the second group.
var noDeNovo = HD.run(Object.assign({}, base, HD.presets.ceiling, { denovo: 0 }));
var fromDeNovo = ceil.b.car[50] - noDeNovo.b.car[50];
check("new mutations are a material part of the residual",
      fromDeNovo / ceil.b.car[50] > 0.25,
      n(fromDeNovo) + " of " + n(ceil.b.car[50]) + " remaining carriers (" +
      (fromDeNovo / ceil.b.car[50] * 100).toFixed(0) + "%) are new expansions");
check("the rest of the residual is the surviving 2026 cohort",
      noDeNovo.b.car[50] < r.a.car[50] * 0.2,
      n(noDeNovo.b.car[50]) + " left of " + n(r.a.car[0]) + " — nobody born after 2027 carries it");

// ---- monotonicity -----------------------------------------------------
var mono = true, prev = Infinity;
[0, 0.25, 0.5, 0.75, 1].forEach(function (c) {
  var v = HD.run(Object.assign({}, base, { cov: c, start: 2030 })).b.car[50];
  if (v > prev + 1e-6) mono = false;
  prev = v;
});
check("more coverage never increases carriers", mono);

// ---- probabilistic sensitivity analysis -------------------------------
var mcOpts = Object.assign({}, base, HD.presets.broad);
var mc = HD.montecarlo(mcOpts, 500, 20260101);
function p1(x) { return (x * 100).toFixed(1) + "%"; }

console.log("\nProbabilistic sensitivity analysis (" + mc.n + " draws, broad eligibility)");
console.log("  2076 carrier reduction    " + p1(mc.reductionCar2076.mid) +
            "   [90% CI " + p1(mc.reductionCar2076.lo) + " – " + p1(mc.reductionCar2076.hi) + "]");
console.log("  2076 diagnosed reduction  " + p1(mc.reductionMan2076.mid) +
            "   [90% CI " + p1(mc.reductionMan2076.lo) + " – " + p1(mc.reductionMan2076.hi) + "]");
console.log("  2076 cases averted        " + n(mc.avertedCases2076.mid) +
            "   [" + n(mc.avertedCases2076.lo) + " – " + n(mc.avertedCases2076.hi) + "]");
console.log("  2076 baseline cases       " + n(mc.bands.aMan.mid[50]) +
            "   [" + n(mc.bands.aMan.lo[50]) + " – " + n(mc.bands.aMan.hi[50]) + "]");
console.log("  top driver                " + mc.drivers[0].param +
            " (rho " + mc.drivers[0].rho.toFixed(2) + ")");
console.log("");

check("the interval brackets its own point estimate",
      mc.reductionCar2076.lo < mc.reductionCar2076.mid &&
      mc.reductionCar2076.mid < mc.reductionCar2076.hi);

check("bands are ordered at every year",
      mc.years.every(function (_, t) {
        return mc.bands.aMan.lo[t] <= mc.bands.aMan.mid[t] &&
               mc.bands.aMan.mid[t] <= mc.bands.aMan.hi[t] &&
               mc.bands.bCar.lo[t] <= mc.bands.bCar.mid[t] &&
               mc.bands.bCar.mid[t] <= mc.bands.bCar.hi[t];
      }));

check("the deterministic run falls inside the band",
      broad.a.man[50] >= mc.bands.aMan.lo[50] && broad.a.man[50] <= mc.bands.aMan.hi[50],
      "point " + n(broad.a.man[50]) + " within " +
      n(mc.bands.aMan.lo[50]) + "–" + n(mc.bands.aMan.hi[50]));

// This is the argument the PSA exists to make: the absolute level is
// uncertain, the proportional policy effect is much less so, because both
// scenarios share the same draws and most of the uncertainty cancels.
var levelSpread = (mc.bands.aMan.hi[50] - mc.bands.aMan.lo[50]) / mc.bands.aMan.mid[50];
var effectSpread = (mc.reductionCar2076.hi - mc.reductionCar2076.lo) / mc.reductionCar2076.mid;
check("the effect ratio is better determined than the level",
      effectSpread < levelSpread,
      "level +/-" + p1(levelSpread / 2) + " vs effect +/-" + p1(effectSpread / 2));

// The strongest honest form of the claim: not an average, a floor over the
// whole sample. If this ever fails, the headline sentence must be softened.
check("every single draw shows a substantial reduction",
      mc.reductionCar2076.min > 0.20,
      "worst of " + mc.n + " draws = " + p1(mc.reductionCar2076.min) +
      ", best = " + p1(mc.reductionCar2076.max));

check("efficacy is the dominant driver", mc.drivers[0].param === "eff",
      "rho = " + mc.drivers[0].rho.toFixed(2));
check("more new mutations means less reducible disease",
      mc.drivers.filter(function (d) { return d.param === "denovo"; })[0].rho < 0,
      "denovo rho = " + mc.drivers.filter(function (d) { return d.param === "denovo"; })[0].rho.toFixed(2));

// A result that moves when the seed changes is a result that needed more
// draws. 500 should hold the median to well inside a percentage point.
var mc2 = HD.montecarlo(mcOpts, 500, 777001);
var drift = Math.abs(mc2.reductionCar2076.mid - mc.reductionCar2076.mid);
check("the estimate is stable across seeds", drift < 0.01,
      "median moves " + p1(drift) + " on a different seed");

check("zero coverage yields no modelled reduction",
      HD.montecarlo({ cov: 0, start: 2035 }, 120, 4242).reductionCar2076.hi < 1e-9);

// ---- narrow exception, the team's actual proposal ---------------------
var narrow = HD.run(Object.assign({}, base, HD.presets.narrow));
console.log("\nNarrow exception (homozygous-parent cases only, ~1% coverage)");
console.log("  2076 cases averted     " + (narrow.a.man[50] - narrow.b.man[50]).toFixed(1));
console.log("  2076 carriers averted  " + (narrow.a.car[50] - narrow.b.car[50]).toFixed(1));
console.log("  (small by construction — the amendment is a necessity argument,");
console.log("   not a population-health program)");

console.log("\n" + "-".repeat(58));
console.log(fails === 0 ? "ALL CHECKS PASS\n" : fails + " CHECK(S) FAILED\n");
process.exit(fails === 0 ? 0 : 1);
