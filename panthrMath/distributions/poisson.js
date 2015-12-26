(function(define) {'use strict';
define(function(require) {

   // No input validation provided.

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the Poisson distribution.
    *
    * @module distributions.poisson
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var lpoisson, lgamma, pgamma, pWrap, qnorm, discInvCdf, inverseCDF;

   lpoisson = require('../basicFunc/lpoisson').lpoisson;
   lgamma = require('../basicFunc/lgamma').lgamma;
   pgamma = require('./gamma').pgamma;
   pWrap = require('../utils').pWrap;
   discInvCdf = require('../utils').discInvCdf;
   qnorm = require('./normal').qnorm;
   inverseCDF = require('../rgen/inverseCDF');

   // density / pdf
   /**
    * TODO
    * @memberof poisson
    */
   function dpois(lambda, logp) {
      logp = logp === true;
      return function(x) {
         return logp ? lpoisson(lambda)(x) : Math.exp(lpoisson(lambda)(x));
      };
   }

   /**
    * TODO
    * @memberof poisson
    */
   function ppois(lambda, lowerTail, logp) {
      lowerTail = lowerTail !== false;
      logp = logp === true;

      if (!(lambda >= 0)) { return function(x) { return NaN; }; }

      return function(x) {
         var ret;

         if (x >= 0 && lambda > 0) {
            return pgamma(Math.floor(x + 1e-10) + 1, 1, !lowerTail, logp)(lambda);
         }

         ret = lowerTail ? 1 : 0;
         if (x < 0) { ret = lowerTail ? 0 : 1; }
         return logp ? Math.log(ret) : ret;
      };
   }

   /**
    * TODO
    * @memberof poisson
    */
   function qpois(lambda, lowerTail, logp) {

      // TODO:  make a separate calculation for !lowerTail

      var goodParams, mu, sigma, gamma;

      logp = logp === true;
      lowerTail = lowerTail !== false;
      goodParams = lambda >= 0 && lambda < Infinity;

      if (!goodParams) { return function(prob) { return NaN; }; }
      if (lambda === 0) {
         return pWrap(lowerTail, logp, function(prob) { return 0; });
      }

      mu = lambda;
      sigma = Math.sqrt(lambda);
      gamma = 1 / sigma;

      return pWrap(true, logp, function(prob) {
         var z, ret;
         z = qnorm(0, 1, lowerTail)(prob); // initial value
         if (z < -10) { z = -10; }
         if (z > 10) { z = 10; }
         ret = Math.floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);
         if (prob === 1) { return lowerTail ? Infinity : 0; }
         if (lowerTail) {
            return discInvCdf(0, 1e10, ret, prob, ppois(lambda));
         }
         return -discInvCdf(-1e10, 0, -ret, prob, function(x) {
            return ppois(lambda, lowerTail)(-x);
         });
      });
   }

   // Using inverseCDF
   /**
    * TODO
    * @memberof poisson
    */
   function rpois(lambda) {
      var mode;

      mode = { val: Math.floor(lambda) };
      mode.prob = Math.exp(-lambda + mode.val * Math.log(lambda) -
         lgamma(mode.val + 1));

      return inverseCDF(
         function getMode() { return mode; },
         function updateLeft() {
            if (this.val === 0) { return false; }
            this.prob *= this.val / lambda;
            this.val -= 1;
            return true;
         },
         function updateRight() {
            this.val += 1;
            this.prob *= lambda / this.val;
            return true;
         }
      );
   }


   return {
      /**
       * Returns an object representing a Poisson distribution, with properties `d`, `p`, `q`, `r`.
       * ```
       * pois(a, b).d(x, logp)            // same as dpois(a, b, logp)(x)
       * pois(a, b).q(x, lowerTail, logp) // same as qpois(a, b, lowerTail, logp)(x)
       * pois(a, b).p(x, lowerTail, logp) // same as ppois(a, b, lowerTail, logp)(x)
       * pois(a, b).r(n)                  // same as rpois(a, b)(n)
       * ```
       * @memberof poisson
       */
      pois: function(lambda) {
         return {
            d: function(x, logp) { return dpois(lambda, logp)(x); },
            p: function(q, lowerTail, logp) {
               return ppois(lambda, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qpois(lambda, lowerTail, logp)(p);
            },
            r: function(n) { return rpois(lambda)(n); }
         };
      },
      dpois: dpois,
      ppois: ppois,
      qpois: qpois,
      rpois: rpois
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
