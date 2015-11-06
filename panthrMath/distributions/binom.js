(function(define) {'use strict';
define(function(require) {

   var C, stirlerr, bd0, pbeta, qnorm, pWrap, discInvCdf;

   C = require('../constants');
   stirlerr = require('../basicFunc/stirlerr').stirlerr;
   bd0 = require('../basicFunc/bd0').bd0;
   pbeta = require('./beta').pbeta;
   qnorm = require('./normal').qnorm;
   pWrap = require('../utils').pWrap;
   discInvCdf = require('../utils').discInvCdf;

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

   function qbinom(n, p, lowerTail, logp) {
      var goodParams, mu, sigma, gamma;

      logp = logp === true;
      lowerTail = lowerTail !== false;
      goodParams = n >= 0 && n < Infinity && n === Math.floor(n)
            && p >= 0 && p <= 1;

      if (!goodParams) { return function(prob) { return NaN; }; }
      if (p === 0) {
         return pWrap(lowerTail, logp,
            function(prob) { return prob === 1 ? n : 0; }
         );
      }
      if (p === 1) {
         return pWrap(lowerTail, logp,
            function(prob) { return prob === 0 ? 0 : n; }
         );
      }
      // 0 < p < 1
      mu = n * p;
      sigma = Math.sqrt(mu * (1 - p));
      gamma = (1 - 2 * p) / sigma;

      return pWrap(lowerTail, logp, function(prob) {
         var z, ret;
         if (prob === 0) { return 0; }
         if (prob === 1) { return n; }
         z = qnorm(0, 1)(prob); // initial value
         ret = Math.floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);
         if (ret > n) { ret = n; }
         return discInvCdf(0, n, ret, prob, pbinom(n, p));
      });
   }

   return {
      binom: function(n, p) {
         return {
            d: function(x, logp) { return dbinom(n, p, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pbinom(n, p, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qbinom(n, p, lowerTail, logp)(p);
            },
            r: function(n) { return rbinom(n, p)(n); }
         };
      },
      dbinom: dbinom,
      pbinom: pbinom,
      qbinom: qbinom
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
