(function(define) {'use strict';
define(function(require) {

   var C, lgamma, lfactorial, lchoose, stirlerr;

   C = require('./constants');
   lgamma = require('./basicFunc/lgamma').lgamma;
   stirlerr = require('./basicFunc/stirlerr').stirlerr;

   // Expects x >= 0;
   lfactorial = function lfactorial(x) {
      return x === 0 ? 0 : lgamma(x + 1);
   };

   lchoose = function lchoose(n, k) {
      return lfactorial(n) - lfactorial(k) - lfactorial(n - k);
   };

   function bd0(x, np) {
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
         throw new Error('bd0(x, np) made it to 10000 ...');
      }
      return x * Math.log(x / np) + np - x;
   }

   // returns the log of the binomial probability
   // Note: the arguments are re-arranged:  lbinomProb(n, p, x)
   // calculates log(binom_prob(x; n, p)) where n > 0 and
   // 0 <= p <= 1.
   // Based on:  "Fast and Accurate Computation of Binomial Probabilities",
   // Loader (2000)
   function lbinomProb(n, p, x) {
      if (p === 0) { return x === 0 ? 0 : -Infinity; }
      if (p === 1) { return x === n ? 0 : -Infinity; }
      if (x === 0) { return n * Math.log(1 - p); }
      if (x === n) { return n * Math.log(p); }
      return stirlerr(n) - stirlerr(x) - stirlerr(n - x) -
         bd0(x, n * p) - bd0(n - x, n * (1 - p)) +
         0.5 * Math.log(n / (C.twopi * x * (n - x)));
   }

   return {
      lfactorial: lfactorial,
      lchoose: lchoose,
      lbinomProb: lbinomProb,
      bd0: bd0
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
