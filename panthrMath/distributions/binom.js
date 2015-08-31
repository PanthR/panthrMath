(function(define) {'use strict';
define(function(require) {

   var C, twopi, Rational, pbinom, qbinom, stirlerr, bd0;

   C = require('../constants');
   twopi = C.twopi;
   stirlerr = require('../basicFunc/stirlerr').stirlerr;
   bd0 = require('../basicFunc/bd0').bd0;

   // returns the log of the binomial probability
   // Note: the arguments are re-arranged:  lbinomProb(n, p, x)
   // calculates log(binom_prob(x; n, p)) where n > 0 and
   // 0 <= p <= 1.
   // Based on:  "Fast and Accurate Computation of Binomial Probabilities",
   // Loader (2000)
   function dbinomLog(n, p) {
      if (p === 0) { return function(x) { return x === 0 ? 0 : -Infinity; }; }
      if (p === 1) { return function(x) { return x === n ? 0 : -Infinity; }; }
      return function(x) {
         if (x === 0) { return n * Math.log(1 - p); }
         if (x === n) { return n * Math.log(p); }
         if (x < 0 || x > n || Math.round(x) !== x) { return -Infinity; }
         return stirlerr(n) - stirlerr(x) - stirlerr(n - x) -
            bd0(x, n * p) - bd0(n - x, n * (1 - p)) +
            0.5 * Math.log(n / (C.twopi * x * (n - x)));
      };
   }

   function dbinom(n, p) {
      var lbinom = dbinomLog(n, p);
      return function(x) {
         return Math.exp(lbinom(x));
      };
   }

   dbinom.log = dbinomLog;

   return {
      dbinom: dbinom,
      pbinom: pbinom,
      qbinom: qbinom
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
