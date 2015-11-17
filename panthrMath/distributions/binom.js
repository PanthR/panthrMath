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
   // Note: the arguments are re-arranged:  lbinomProb(size, p, x)
   // calculates log(binom_prob(x; size, p)) where size > 0 and
   // 0 <= p <= 1.
   // Based on:  "Fast and Accurate Computation of Binomial Probabilities",
   // Loader (2000)
   function dbinomLog(size, p) {
      if (p === 0) { return function(x) { return x === 0 ? 0 : -Infinity; }; }
      if (p === 1) { return function(x) { return x === size ? 0 : -Infinity; }; }
      return function(x) {
         if (x === 0) { return size * Math.log(1 - p); }
         if (x === size) { return size * Math.log(p); }
         if (x < 0 || x > size || Math.round(x) !== x) { return -Infinity; }
         return stirlerr(size) - stirlerr(x) - stirlerr(size - x) -
            bd0(x, size * p) - bd0(size - x, size * (1 - p)) +
            0.5 * Math.log(size / (C.twopi * x * (size - x)));
      };
   }

   function dbinom(size, p, logp) {
      var lbinom;

      logp = logp === true;
      lbinom = dbinomLog(size, p);

      return function(x) {
         return logp ? lbinom(x) : Math.exp(lbinom(x));
      };
   }

   function pbinom(size, p, lowerTail, logp) {
      var goodParams;

      logp = logp === true;
      lowerTail = lowerTail !== false;

      goodParams = size >= 0 && size < Infinity && size === Math.floor(size)
            && p >= 0 && p <= 1;
      if (!goodParams) { return function(x) { return NaN; }; }
      return function(x) {
         var res;

         if (isNaN(x)) { return NaN; }
         res = lowerTail === x < 0 ? -Infinity : 0;
         x = Math.floor(x);
         if (x >= 0 && x < size) {
            res = pbeta(x + 1, size - x, !lowerTail, true)(p);
         }

         return logp ? res : Math.exp(res);
      };
   }

   function qbinom(size, p, lowerTail, logp) {

      // TODO:  make a separate calculation for !lowerTail

      var goodParams, mu, sigma, gamma;

      logp = logp === true;
      lowerTail = lowerTail !== false;
      goodParams = size >= 0 && size < Infinity && size === Math.floor(size)
            && p >= 0 && p <= 1;

      if (!goodParams) { return function(prob) { return NaN; }; }
      if (p === 0) {
         return pWrap(lowerTail, logp,
            function(prob) { return prob === 1 ? size : 0; }
         );
      }
      if (p === 1) {
         return pWrap(lowerTail, logp,
            function(prob) { return prob === 0 ? 0 : size; }
         );
      }
      // 0 < p < 1
      mu = size * p;
      sigma = Math.sqrt(mu * (1 - p));
      gamma = (1 - 2 * p) / sigma;

      return pWrap(true, logp, function(prob) {
         var z, ret;
         // if (prob === 0) { return 0; }
         // if (prob === 1) { return size; }
         z = qnorm(0, 1, lowerTail)(prob); // initial value
         ret = Math.floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);
         if (ret > size) { ret = size; }
         if (prob === 1) { return lowerTail ? size : 0; }
         if (lowerTail) {
            return discInvCdf(0, size, ret, prob, pbinom(size, p));
         }
         return -discInvCdf(-size, 0, -ret, prob, function(x) {
            return pbinom(size, p, lowerTail)(-x);
         });
      });
   }

   return {
      binom: function(size, p) {
         return {
            d: function(x, logp) { return dbinom(size, p, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pbinom(size, p, lowerTail, logp)(q);
            },
            q: function(prob, lowerTail, logp) {
               return qbinom(size, p, lowerTail, logp)(prob);
            },
            r: function(n) { return rbinom(size, p)(n); }
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
