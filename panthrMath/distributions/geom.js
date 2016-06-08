(function(define) {
'use strict';
define(function(require) {

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the geometric distribution with parameter the probability $p$, which
    * is defined by the pdf
    * $$ f(x; p) = p (1-p)^x$$
    * where $0 < p \leq 1$ and $x \in\left\\{0,1,2,\ldots\right\\}$.
    *
    * `dgeom` provides access to this probability density function,
    * `pgeom` to the cumulative distribution function, `qgeom` to the
    * quantile function (inverse cdf)
    * and `rgeom` to random deviates.
    *
    * Finally, you can use `geom` to obtain an object
    * representing the geometric distribution for a given value of
    * the probability parameter $p$.
    *
    * @module distributions.geometric
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var log1p, expm1, rexp, rpois, dbinomLog, utils;

   rexp = require('./exp').rexp;
   rpois = require('./poisson').rpois;
   log1p = require('../basicFunc/log1p').log1p;
   expm1 = require('../basicFunc/expm1').expm1;
   dbinomLog = require('../basicFunc/dbinomLog').dbinomLog;
   utils = require('../utils');

   /**
    * Evaluates the geometric distribution's density function at `x`:
    * $$\textrm{dgeom}(p)(x) = p(1-p)^x$$
    * for $x \in\left\\{0,1,2,\ldots\right\\}$.
    *
    * Expects the probability $0 < p\leq 1$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the
    * logarithm of the result.
    *
    * The implementation uses `dbinom`.
    *
    * @fullName dgeom(prob, logp)(x)
    * @memberof geometric
    */
   function dgeom(prob, logp) {
      logp = logp === true;

      if (prob <= 0 || prob > 1 || utils.hasNaN(prob)) {
         return function(x) { return NaN; };
      }

      return function(x) {
         if (utils.hasNaN(x)) { return NaN; }
         if (x < 0 || Math.round(x) !== x) {
            return logp ? -Infinity : 0;
         }
         return logp ? Math.log(prob) + dbinomLog(x, prob)(0)
                     : prob * Math.exp(dbinomLog(x, prob)(0));
      };
   }

   /**
    * Evaluates the lower-tail cdf at `x` for the geometric distribution:
    * $$\textrm{pgeom}(p)(x) = 1 - (1-p)^{x+1}$$
    *
    * Expects the probability $0 < p \leq 1$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability instead.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    *
    * The implementation uses `pbinom`.
    *
    * @fullName pgeom(prob, lowerTail, logp)(x)
    * @memberof geometric
    */
   function pgeom(prob, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (prob <= 0 || prob > 1 || utils.hasNaN(prob)) {
         return function(x) { return NaN; };
      }
      if (prob === 1) {
         return function(x) {
            var ret;

            if (utils.hasNaN(x)) { return NaN; }
            ret = x < 0 ? 0 : 1;
            ret = lowerTail ? ret : 1 - ret;
            return logp ? Math.log(ret) : ret;
         };
      }

      return function(x) {
         var ret;

         if (utils.hasNaN(x)) { return NaN; }

         x = Math.floor(x + 1e-14);

         if (x < 0) {
            ret = lowerTail ? 0 : 1;
            return logp ? Math.log(ret) : ret;
         }

         // basic behavior
         ret = log1p(-prob) * (x + 1); // log ((1-p)^(x+1));
         if (lowerTail) {
            return logp ? log1p(-Math.exp(ret)) : -expm1(ret);
         }
         // else right tail
         return logp ? ret : Math.exp(ret);
      };
   }

   /**
    * Evaluates the geometric distribution's quantile function
    * (inverse cdf) at `p`:
    * $$\textrm{qexp}(\textrm{prob})(p) = x \textrm{ such that } \textrm{prob}(X \leq x) = p$$
    * where $X$ is a random variable with the geometric distribution.
    *
    * Expects the probability parameter $0 < \textrm{prob} \leq 1$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
    * interpreted as an upper tail probability (returns
    * $x$ such that $\textrm{prob}(X > x) = p)$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
    * the logarithm of the desired probability.
    * @fullName qgeom(prob, lowerTail, logp)(p)
    * @memberof geometric
    */
   function qgeom(prob, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (prob <= 0 || prob > 1) { return function(x) { return NaN; }; }

      return function(p) {
         var res;

         if (logp) {
            if (p > 0) { return NaN; }
            res = lowerTail ? Math.log(-expm1(p)) : p;
         } else {
            // else, not log
            if (!(p >= 0 && p <= 1)) { return NaN; }
            res = lowerTail ? log1p(-p) : Math.log(p);
         }
         res = res / log1p(-prob) - 1;

         return Math.max(0, Math.ceil(res - 1e-12));
      };
   }

   /**
    * Returns a random variate from the geometric distribution.
    *
    * Expects the probability parameter $0 < p \leq 1$.
    *
    * Following R's code (rgeom.c)
    * @fullName rgeom(prob)()
    * @memberof geometric
    */
   function rgeom(prob) {
      var scale;

      if (prob <= 0 || prob > 1) { return function() { return NaN; }; }

      scale = prob / (1 - prob);

      return function() {
         return rpois(rexp(scale)())();
      };
   }

   return {
      /**
       * Returns an object representing an geometric distribution with a
       * given probability parameter $0 < p \leq 1$,
       * with properties `d`, `p`, `q`, `r`.
       * ```
       * geom(prob).d(x, logp)            // same as dgeom(prob, logp)(x)
       * geom(prob).p(x, lowerTail, logp) // same as pgeom(prob, lowerTail, logp)(x)
       * geom(prob).q(x, lowerTail, logp) // same as qgeom(prob, lowerTail, logp)(x)
       * geom(prob).r()                   // same as rgeom(prob)()
       * ```
       * @memberof geometric
       */
      geom: function(prob) {
         return {
            d: function(x, logp) { return dgeom(prob, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pgeom(prob, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qgeom(prob, lowerTail, logp)(p);
            },
            r: function() { return rgeom(prob)(); }
         };
      },
      dgeom: dgeom,
      pgeom: pgeom,
      qgeom: qgeom,
      rgeom: rgeom
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
