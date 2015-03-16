(function(define) {'use strict';
define(function(require) {

   var C, lgamma, lfactorial, lchoose, lbinomProb, stirlerr;

   C = require('./constants');

   lgamma = (function() {
      // z Must be positive
      var cs = [
         0.99999999999980993227684700473478,
         676.520368121885098567009190444019, -1259.13921672240287047156078755283,
         771.3234287776530788486528258894, -176.61502916214059906584551354,
         12.507343278686904814458936853, -0.13857109526572011689554707,
         9.984369578019570859563e-6, 1.50563273514931155834e-7
      ];
      return function(z) {
         var t, i, ser;

         t = z + 7.5;
         t -= (z + 0.5) * Math.log(t);
         ser = cs[0];
         for (i = 1; i < cs.length; i += 1) {
            ser += cs[i] / (z + i);
         }

         return -t + Math.log(C.sqrt2pi * ser / z);
      };
   }());

   // Expects x >= 0;
   lfactorial = function lfactorial(x) {
      return x === 0 ? 0 : lgamma(x + 1);
   };

   lchoose = function lchoose(n, k) {
      return lfactorial(n) - lfactorial(k) - lfactorial(n - k);
   };

   // error term in Stirling's approximation
   // log(n!) - log( sqrt(2*pi*n) * (n/e)^n )
   stirlerr = (function() {
      var cs, precomputed, N;
      // coefficients in Stirling's expansion for log(Gamma)
      cs = [1 / 12, 1 / 360, 1 / 1260, 1 / 1680, 1 / 1188];
      precomputed = [NaN];
      function compute(n) {
         var nsq = n * n;
         return (cs[0] - (cs[1] - (cs[2] - (cs[3] - cs[4] / nsq) / nsq) / nsq ) / nsq) / n;
      }
      for (N = 1; N < 16; N += 1) {
         precomputed[N] = compute(N);
      }
      return function(n) { return precomputed[n] || compute(n); };
   }());

   // returns the log of the binomial probability
   // Note: the arguments are re-arranged:  lbinomProb(n, p, x)
   // calculates log(binom_prob(x; n, p)) where n > 0 and
   // 0 <= p <= 1.
   // Based on:  "Fast and Accurate Computation of Binomial Probabilities",
   // Loader (2000)
   lbinomProb = (function() {
      function dev(np, x) {
         var ej, j, s, sTemp, v;
         if (Math.abs(x - np) < 0.1 * (x + np)) {
            v = (x - np) / (x + np);
            s = (x - np) * v;
            ej = 2 * x * v;
            for (j = 1; j < 10000; j += 1) {
               ej *= v * v;
               sTemp = s + ej / (2 * j + 1);
               if (sTemp === s) {
                  return sTemp;
               }
               s = sTemp;
            }
            throw new Error('dev(np, x) made it to 10000 ...');
         }
         return x * Math.log(x / np) + np - x;
      }

      return function(n, p, x) {
         if (p === 0) { return x === 0 ? 0 : -Infinity; }
         if (p === 1) { return x === n ? 0 : -Infinity; }
         if (x === 0) { return n * Math.log(1 - p); }
         if (x === n) { return n * Math.log(p); }
         return stirlerr(n) - stirlerr(x) - stirlerr(n - x) -
            dev(n * p, x) - dev(n * (1 - p), n - x) +
            0.5 * Math.log(n / (C.twopi * x * (n - x)));
      };
   }());


   return {
      lgamma: lgamma,
      lfactorial: lfactorial,
      lchoose: lchoose,
      lbinomProb: lbinomProb,
      stirlerr: stirlerr
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
