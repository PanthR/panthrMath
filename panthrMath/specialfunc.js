(function(define) {'use strict';
define(function(require) {

   var C, lgamma, lfactorial, lchoose, stirlerr, bd0;

   C = require('./constants');
   lgamma = require('./basicFunc/lgamma').lgamma;
   stirlerr = require('./basicFunc/stirlerr').stirlerr;
   bd0 = require('./basicFunc/bd0').bd0;

   // Expects x >= 0;
   lfactorial = function lfactorial(x) {
      return x === 0 ? 0 : lgamma(x + 1);
   };

   lchoose = function lchoose(n, k) {
      return lfactorial(n) - lfactorial(k) - lfactorial(n - k);
   };

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
      lbinomProb: lbinomProb
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
