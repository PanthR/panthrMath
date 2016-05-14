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
var nks = 1500;
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
   }

   // Return proportion of small p-values. On a properly calibrated r-distr
   // this should be 0-20%.
   return pvalues.reduce(function(acc, v) {
         return v < 0.05 ? (acc + 1) : acc;
      }, 0) / reps;
}

/*
 * This will do a Kolmogorov-Smirnov goodness-of-fit test.
 */
function testCont(distrObj) {
   var pvalues, rep, observed, ksdiffs, dstat, i;

   pvalues = [];

   for (rep = 0; rep < reps; rep += 1) {
      observed = [];
      for (i = 0; i < nks; i += 1) {
         observed.push(distrObj.r());
      }
      observed.sort(function(x, y) { return x - y; });
      ksdiffs = observed.map(function(x, i) {
         return distrObj.p(x) - i / nks;
      });
      dstat = ksdiffs.reduce(function(acc, x) {
         return Math.max(acc, Math.abs(x), Math.abs(1 / nks - x));
      }, -1);

      pvalues.push(ks(nks, dstat));
   }

   return pvalues.reduce(function(acc, v) {
      return v < 0.05 ? (acc + 1) : acc;
   }, 0) / reps;

}


// BEGIN KS IMPLEMENTATION

 /* Compute Kolmogorov's distribution.
    Code published in
 George Marsaglia and Wai Wan Tsang and Jingbo Wang (2003),
 "Evaluating Kolmogorov's distribution".
 Journal of Statistical Software, Volume 8, 2003, Issue 18.
 URL: http://www.jstatsoft.org/v08/i18/.

 We are returning the upper-tail p-value
 */
function ks(n, d)
{
   var k, m, i, j, g, h, s, H, Q, basicArr;

   /*
      The faster right-tail approximation is omitted here.
   */
   s = d * d * n;
   if (s > 7.24 || (s > 3.76 && n > 99)) {
      return 2 * Math.exp(-(2.000071+.331/Math.sqrt(n)+1.409/n)*s);
   }
   k = Math.floor(n * d) + 1;
   m = 2 * k - 1;
   h = k - n * d;
   basicArr = seq(1, m * m);
   H = basicArr.map(function() { return 0; });
   Q = basicArr.map(function() { return 0; });
   for (i = 0; i < m; i++) {
       for (j = 0; j < m; j++) {
         H[i * m + j] = i - j + 1 < 0 ? 0 : 1;
       }
   }
   // Fixing "numerators"
   for (i = 0; i < m; i++) {
       H[i * m] -= Math.pow(h, i + 1);
       H[(m - 1) * m + i] -= Math.pow(h, m - i);
   }
   // Adjust the bottom left corner if h > 1/2
   H[(m - 1) * m] += 2 * h - 1 > 0 ? Math.pow(2 * h - 1, m) : 0;
   // Fix "denominators"
   // TODO: Maybe can do better?
   for (i = 0; i < m; i++) {
       for (j = 0; j < i + 1; j++) {
            for (g = 1; g <= i - j + 1; g++) {
               H[i * m + j] /= g;
            }
       }
   }

   // Q will be an object with "arr" and "exp" properties
   Q = m_power({ arr: H, exp: 0 }, m, n);

   s = Q.arr[(k - 1) * m + k - 1];
   for(i = 1; i <= n; i++) {
       s = s * i / n;
       if (s < 1e-140) {
         s *= 1e140;
         Q.exp -= 140;
       }
   }
   return 1 - s * Math.pow(10, Q.exp);
}

// Multiplication of two mxm matrices.
function m_multiply(A, B, m) {
   var i, j, k, C;
   C = [];
   for (i = 0; i < m; i++) {
      for (j = 0; j < m; j++) {
         C[i * m + j] = 0;
         for(k = 0; k < m; k++) {
            C[i * m + j] += A[i * m + k] * B[k * m + j];
         }
      }
   }
   return C;
}

// Computes the "n"-th power of an mxm matrix
// arrObj is an object with properties arr and exp
// Returns an object with the same properties.
function m_power(arrObj, m, n) {

   var halfObj, retObj, i;

   if(n == 1) {
      return arrObj;
   }

   halfObj = m_power(arrObj, m, Math.floor(n / 2));

   retObj = { arr: m_multiply(halfObj.arr, halfObj.arr, m),
              exp: 2 * halfObj.exp };

   if (n % 2 == 1) {
      retObj = { arr: m_multiply(arrObj.arr, retObj.arr, m),
                 exp: arrObj.exp + retObj.exp };
    }

    if (retObj.arr[(m / 2) * m + (m / 2)] > 1e140) {
      for(i = 0; i < m * m; i++) {
         retObj.arr[i] *= 1e-140;
      }
      retObj.exp += 140;
    }

    return retObj;
}

// END KS IMPLEMENTATION


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
   describe("for negative binomial", function() {
      seq(1, 10).forEach(function() {
         it("repeated test", function() {
            var n, p, pvalues, xs, obj, upper, lower;
            n = 5 + Math.floor(Math.random() * 10);
            p = main.runif(0.5, 0.8)();
            obj = main.nbinom(n, p);
            upper = n;
            while ( obj.d(upper) >= 1e-4) { upper += 1; }

            expect(testDiscrete(obj, seq(0, upper), 2)).to.be.below(0.2);
            // Providing a "closeby" distribution
            obj.r = main.nbinom(n, p * 0.9).r;
            expect(testDiscrete(obj, seq(0, n), 2)).to.be.above(0.2);
         });
      });
   });
   describe("for geometric", function() {
      seq(1, 20).forEach(function() {
         it("repeated test", function() {
            var prob, pvalues, xs, obj, n, p;
            prob = main.runif(0.1, 0.9)();
            // Finding the point where right tail small enough
            n = 0;
            p = prob;
            while (p >= 1e-4) { n += 1; p *= 1 - prob; }
            expect(testDiscrete(main.geom(prob), seq(0, n), 1)).to.be.below(0.2);
            // Providing a "closeby" distribution
            obj = main.geom(prob);
            obj.r = main.geom(prob * 1.05).r;
            expect(testDiscrete(obj, seq(0, n), 1)).to.be.above(0.2);
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
   describe("for normal", function() {
      it("repeated test", function() {
         var mu = Math.random() * 5;
         var sigma = 1 + Math.random() * 3;
         var distr = main.normal(mu, sigma);
         expect(testCont(distr)).to.be.below(0.2);
         distr.r = main.normal(mu, 1.1 * sigma).r;
         expect(testCont(distr)).to.be.above(0.2);
         distr.r = main.normal(mu + 0.1 * sigma, sigma).r;
         expect(testCont(distr)).to.be.above(0.2);
      });
   });
   describe("for lognormal", function() {
      seq(1, 5).forEach(function() {
         it("repeated test", function() {
            var mu = Math.random() * 5;
            var sigma = 1 + Math.random() * 3;
            var distr = main.lognormal(mu, sigma);
            expect(testCont(distr)).to.be.below(0.2);
            distr.r = main.lognormal(mu, 1.1 * sigma).r;
            expect(testCont(distr)).to.be.above(0.2);
            distr.r = main.lognormal(mu + 0.1 * sigma, sigma).r;
            expect(testCont(distr)).to.be.above(0.2);
         });
      });
   });
   describe("for exponential", function() {
      it("repeated test", function() {
         var lambda = 0.1 + Math.random() * 10;
         var distr = main.expdistr(lambda);
         expect(testCont(distr)).to.be.below(0.2);
         distr.r = main.expdistr(1.05 * lambda).r;
         expect(testCont(distr)).to.be.above(0.2);
      });
   });
   describe("for uniform", function() {
      it("repeated test", function() {
         var min = Math.random() * 10;
         var max = min + Math.random() * 20;
         var distr = main.unif(min, max);
         expect(testCont(distr)).to.be.below(0.2);
         distr.r = main.unif(min, 1.1 * max).r;
         expect(testCont(distr)).to.be.above(0.2);
         distr.r = main.unif(1.1 * min, max).r;
         expect(testCont(distr)).to.be.above(0.2);
      });
   });
   describe("for t", function() {
      seq(1, 10).forEach(function() {
         it("repeated test", function() {
            var df = Math.random() * 30;
            var distr = main.tdistr(df);
            expect(testCont(distr)).to.be.below(0.2);
         });
      });
   });
   describe("for F", function() {
      seq(1, 10).forEach(function() {
         it("repeated test", function() {
            var df1 = Math.random() * 30;
            var df2 = Math.random() * 30;
            var distr = main.fdistr(df1, df2);
            expect(testCont(distr)).to.be.below(0.2);
         });
      });
   });
   describe("for gamma", function() {
      it("repeated test", function() {
         var s = 1 + Math.random() * 4;
         var a = 1 + Math.random() * 4;
         var distr = main.gammadistr(a, s);
         expect(testCont(distr)).to.be.below(0.2);
         distr.r = main.gammadistr(a, 1.1 * s).r;
         expect(testCont(distr)).to.be.above(0.2);
         distr.r = main.gammadistr(1.1 * a, s).r;
         expect(testCont(distr)).to.be.above(0.2);
      });
   });
   describe("for beta", function() {
      it("a param > 1 but close to 1", function() {
         var a = 1 + Math.random() * 0.1;
         var b = 2 + Math.random() * 6;
         var distr = main.betadistr(a, b);
         expect(testCont(distr)).to.be.below(0.2);
         distr.r = main.betadistr(a, 1.1 * b).r;
         expect(testCont(distr)).to.be.above(0.2);
         distr.r = main.betadistr(1.1 * a, b).r;
         expect(testCont(distr)).to.be.above(0.2);
      });
      it("both params > 1", function() {
         var a = 2 + Math.random() * 6;
         var b = 2 + Math.random() * 6;
         var distr = main.betadistr(a, b);
         expect(testCont(distr)).to.be.below(0.2);
         distr.r = main.betadistr(a, 1.1 * b).r;
         expect(testCont(distr)).to.be.above(0.2);
         distr.r = main.betadistr(1.1 * a, b).r;
         expect(testCont(distr)).to.be.above(0.2);
      });
      it("one param < 1", function() {
         var a = Math.random();
         var b = 1 + Math.random() * 6;
         distr = main.betadistr(a, b);
         expect(testCont(distr)).to.be.below(0.2);
      });
      it("both params < 1", function() {
         // Values between (0,1) have different algorithm
         var a = 0.5 + Math.random() * 0.5;
         var b = 0.5 + Math.random() * 0.5;
         var distr = main.betadistr(a, b);
         expect(testCont(distr)).to.be.below(0.2);
      });
      it("one param close to 0", function() {
         // Values between (0,1) have different algorithm
         var a = Math.random() * 0.1;
         var b = Math.random() * 0.5 + 0.5;
         var distr = main.betadistr(a, b);
         expect(testCont(distr)).to.be.below(0.2);
      });
   });
   describe("for chi-squared", function() {
      it("repeated test", function() {
         var df, distr;
         for (var i = 0; i < 5; i += 1) {
            df = Math.random() * 60;
            distr = main.chisq(df);
            expect(testCont(distr)).to.be.below(0.2);
         }
      });
   });
});
