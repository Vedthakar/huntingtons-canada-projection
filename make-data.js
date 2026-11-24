#!/usr/bin/env node
/* Regenerates the published CSVs.  Run:  node make-data.js
   Everything under data/ is OUTPUT of the model plus the parameter sheet that
   goes into it. There is no input dataset: nothing here was fitted. */

var fs = require("fs");
var HD = require("./model.js");
var base = HD.defaults;

if (!fs.existsSync("data")) fs.mkdirSync("data");

function csv(rows) {
  return rows.map(function (r) {
    return r.map(function (c) {
      c = c === null || c === undefined ? "" : String(c);
      return /[",\n]/.test(c) ? '"' + c.replace(/"/g, '""') + '"' : c;
    }).join(",");
  }).join("\n") + "\n";
}
function r2(x) { return Math.round(x * 100) / 100; }

// ---- 1. the parameter sheet -------------------------------------------
var params = [
  ["parameter", "symbol", "value", "units", "prior_distribution", "prior_range", "source"],
  ["Diagnosed prevalence, 2026", "", 13.7, "per 100,000", "normal (sd 1.3)", "9 - 19",
   "Fisher & Hayden, British Columbia population study, Movement Disorders, 2014"],
  ["Mendelian transmission", "", 0.5, "probability per conception", "fixed", "",
   "Autosomal dominant inheritance; not uncertain"],
  ["De novo share of incidence", "phi", 0.06, "share", "uniform", "0.04 - 0.10",
   "Standard 5-10% range for new expansions from intermediate alleles (27-35 CAG)"],
  ["Per-embryo editing efficacy", "epsilon", 0.90, "share", "triangular", "0.70 - 0.99",
   "Scenario judgement: pessimistic to near-perfect editing"],
  ["Full-penetrance share (>=40 CAG)", "pi", 0.95, "share", "uniform", "0.92 - 0.98",
   "Caron, Wright & Hayden, GeneReviews, 2020; Kay et al., Neurology 87(3), 2016"],
  ["Onset distribution median", "", 45, "years", "triangular", "40 - 50",
   "Caron, Wright & Hayden, GeneReviews, 2020"],
  ["Onset distribution scale", "", 8, "years", "uniform", "6.5 - 9.5",
   "Width of the logistic onset curve"],
  ["Excess hazard once symptomatic", "", 0.055, "per year", "uniform", "0.045 - 0.067",
   "Implies roughly 15-22 year median survival from onset"],
  ["Background mortality", "m(a)", "0.00022*exp(0.085a)+0.0005", "per year", "fixed", "",
   "Gompertz fit to Canadian all-cause mortality"],
  ["Total fertility rate", "TFR", 1.26, "births per woman", "triangular", "1.05 - 1.55",
   "Statistics Canada, 2023"],
  ["Population, 2026", "", 41500000, "people", "fixed", "",
   "Statistics Canada population estimate"],
  ["Annual population growth", "", 0.006, "per year", "triangular", "0.002 - 0.011",
   "Statistics Canada projections, low to high scenarios"],
  ["Homozygosity among HD cases", "", 0.003, "share", "not sampled", "",
   "Sengupta & Dasgupta, Clinical Case Reports, 2026; sets the narrow preset's ~1% coverage"]
];
fs.writeFileSync("data/parameters.csv", csv(params));

// ---- 2. year-by-year projections for every preset ----------------------
var rows = [["preset", "coverage", "policy_start_year", "year", "population",
             "scenario", "diagnosed", "carriers", "diagnosed_per_100k",
             "carriers_per_100k", "births_prevented_cumulative"]];
Object.keys(HD.presets).forEach(function (name) {
  var p = HD.presets[name];
  var out = HD.run(Object.assign({}, base, p));
  out.years.forEach(function (y, t) {
    [["prohibition retained", "a"], ["AHRA amended", "b"]].forEach(function (s) {
      var d = out[s[1]];
      rows.push([name, p.cov, p.start, y, Math.round(out.pop[t]), s[0],
                 r2(d.man[t]), r2(d.car[t]),
                 r2(d.man[t] / out.pop[t] * 1e5), r2(d.car[t] / out.pop[t] * 1e5),
                 s[1] === "b" ? r2(out.prevented[t]) : 0]);
    });
  });
});
fs.writeFileSync("data/projections.csv", csv(rows));

// ---- 3. Monte Carlo bands ---------------------------------------------
var mc = HD.montecarlo(Object.assign({}, base, HD.presets.broad), 500, 20260101);
var bandRows = [["year", "scenario", "measure", "p05", "p50", "p95"]];
[["aMan", "prohibition retained", "diagnosed"],
 ["bMan", "AHRA amended", "diagnosed"],
 ["aCar", "prohibition retained", "carriers"],
 ["bCar", "AHRA amended", "carriers"]].forEach(function (k) {
  var b = mc.bands[k[0]];
  mc.years.forEach(function (y, t) {
    bandRows.push([y, k[1], k[2], r2(b.lo[t]), r2(b.mid[t]), r2(b.hi[t])]);
  });
});
fs.writeFileSync("data/monte-carlo-bands.csv", csv(bandRows));

// ---- 4. driver ranking -------------------------------------------------
var drv = [["input", "spearman_rho_vs_2076_carrier_reduction"]];
mc.drivers.forEach(function (d) { drv.push([d.param, d.rho.toFixed(4)]); });
fs.writeFileSync("data/sensitivity-drivers.csv", csv(drv));

console.log("data/parameters.csv            " + (params.length - 1) + " parameters");
console.log("data/projections.csv           " + (rows.length - 1) + " rows");
console.log("data/monte-carlo-bands.csv     " + (bandRows.length - 1) + " rows (500 draws, seed 20260101)");
console.log("data/sensitivity-drivers.csv   " + (drv.length - 1) + " inputs");
