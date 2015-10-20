(function(define) {'use strict';
define(function(require) {

   var C, qbinom, stirlerr, bd0, pbeta;

   C = require('../constants');
   stirlerr = require('../basicFunc/stirlerr').stirlerr;
   bd0 = require('../basicFunc/bd0').bd0;
   pbeta = require('./beta').pbeta;

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

   function dbinom(n, p, logp) {
      var lbinom;

      logp = logp === true;
      lbinom = dbinomLog(n, p);

      return function(x) {
         return logp ? lbinom(x) : Math.exp(lbinom(x));
      };
   }

   function pbinom(n, p, lowerTail, logp) {
      var goodParams;

      logp = logp === true;
      lowerTail = lowerTail !== false;

      goodParams = n >= 0 && n < Infinity && n === Math.floor(n)
            && p >= 0 && p <= 1;
      if (!goodParams) { return function(x) { return NaN; }; }
      return function(x) {
         var res;

         if (isNaN(x)) { return NaN; }
         res = lowerTail === x < 0 ? -Infinity : 0;
         x = Math.floor(x);
         if (x >= 0 && x < n) {
            res = pbeta(x + 1, n - x, !lowerTail, true)(p);
         }

         return logp ? res : Math.exp(res);
      };
   }

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
