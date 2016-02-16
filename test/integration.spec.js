/* Integration tests to ensure random number generators produce the correct
 * distributions.
 * This likely does not check for "proper" "pseudo-random" behavior.
 *
 * Simply uses chi-squared test to compare the distribution produced by
 * the random variate generators to the theoretical one.
 */
var main = require('..');
var chai = require('chai');
var expect = chai.expect;
var utils = require('../panthrMath/utils');
var chisq = main.chisq;
var precision = 1e-8;

var n = 10000;    // Number of random variates generated for each test
var reps = 100;    // Number of times to repeat test

// Returns the array [a, a+1, ..., b]
function seq(a, b) {
   var ret;
   ret = [];
   while (a <= b) {
      ret.push(a);
      a += 1;
   }
   return ret;
}

// Remove the x values that are have low expected frequencies and also low actual
// requencies
function cleanUp(xs, d) {
   return xs.filter(function(x) {
      return d(x) > 1e-3;
   });
}

function getPvalue(observed, predicted, dfs) {
   var testStat;

   testStat = predicted.reduce(function(acc, p, i) {
      var ratio = ((observed[i] / n) - p) / p;
      return acc + p * ratio * ratio;
   }, 0);
   testStat *= n;
   // console.log("testStat", testStat, "dfs", dfs);
   return chisq(dfs).p(testStat, false);
}
/*
   * Returns the array of `reps` many p-values from the test
   * `distrObj` is the distribution object (with r, p, q, d properties)
   * `xs` is the array of possible values for the distribution.
   * The function will filter these to something closer to the
   * distribution's "support" before applying the chi.sq test.
*/
function testDiscrete(distrObj, xs, nparams) {
   var pvalues, observed, predicted, rep, i, variate, filterFun;

   pvalues = [];
   predicted = xs.map(function(x) { return distrObj.d(x); });

   for (rep = 0; rep < reps; rep += 1) {
      observed = xs.map(function(x) { return 0; });

      for (i = 0; i < n; i += 1) {
         variate = distrObj.r();
         xs.forEach(function(x, i) {
            if (variate == x) { observed[i] += 1; }
         });
      }
      filterFun = function(v, i) {
         return observed[i] > 4 && predicted[i] > 4 / n;
      };
      pvalues.push(getPvalue(observed.filter(filterFun),
                             predicted.filter(filterFun),
                             xs.filter(filterFun).length - nparams - 1));
      // console.log(xs.filter(filterFun));
   }

   // Return proportion of small p-values. On a properly calibrated r-distr
   // this should be 0-20%.
   return pvalues.reduce(function(acc, v) {
         return v < 0.05 ? (acc + 1) : acc;
      }, 0) / reps;
}

describe("(slow) Integration test", function() {
   describe("for binomial", function() {
      seq(1, 20).forEach(function() {
         it("repeated test", function() {
            var n, p, pvalues, xs, obj;
            n = 20 + Math.floor(Math.random() * 30);
            p = main.runif(0.1, 0.9)();
            expect(testDiscrete(main.binom(n, p), seq(0, n), 2)).to.be.below(0.2);
            // Providing a "closeby" distribution
            obj = main.binom(n, p);
            obj.r = main.binom(n, p * 1.01).r;
            expect(testDiscrete(obj, seq(0, n), 2)).to.be.above(0.2);
         });
      });
   });
   describe("for Poisson", function() {
      seq(1, 20).forEach(function() {
         it("repeated test", function() {
            var lambda, range, pvalues, xs, obj;
            lambda = main.runif(0.2, 4)();
            range = seq(0, lambda + 10 * Math.sqrt(lambda));
            expect(testDiscrete(main.pois(lambda), range, 1)).to.be.below(0.2);
            // Providing a "closeby" distribution
            obj = main.pois(lambda);
            obj.r = main.pois(lambda * 1.02).r;
            expect(testDiscrete(obj, range, 1)).to.be.above(0.2);
         });
      });
   });
   describe("for finite", function() {
      it("repeated test", function() {
         var setObj, dObj, pvalues, xs, obj;
         setObj = {
            xs: [1, 5, 7, 10],
            ws: [0.1, 0.2, 0.3, 0.4]
         };
         dObj = main.finite(setObj);
         expect(testDiscrete(dObj, setObj.xs, 0)).to.be.below(0.2);
         // Providing a "closeby" distribution
         setObj.ws[3] = 0.43;
         dObj.r = main.finite(setObj).r;
         expect(testDiscrete(dObj, setObj.xs, 0)).to.be.above(0.2);
      });
   });
});
