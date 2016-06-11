(function(define) {
'use strict';
define(function(require) {

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the continuous exponential distribution with *rate* $\lambda$, which
    * is defined by the pdf
    * $$ f(x; \lambda) = \lambda e^{-\lambda x}$$
    * where $\lambda > 0$ and $x \in [0, \infty)$.
    *
    * `dexp` provides access to this probability density function,
    * `pexp` to the cumulative distribution function, `qexp` to the
    * quantile function (inverse cdf)
    * and `rexp` to random deviates.
    *
    * Finally, you can use `expdistr` to obtain an object
    * representing the exponential distribution for a given value of
    * the *rate* $\lambda$.
    *
    * @module distributions.exponential
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var sexp, log1p, expm1, rand, utils;

   log1p = require('../basicFunc/log1p').log1p;
   expm1 = require('../basicFunc/expm1').expm1;
   rand = require('../rgen/rgen').random;
   utils = require('../utils');

   /**
    * Evaluates the exponential distribution's density function at `x`:
    * $$\textrm{dexp}(\lambda)(x) = \lambda e^{-\lambda x}$$
    * for $x \geq 0$.
    *
    * Expects the rate $\lambda > 0$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the
    * logarithm of the result.
    *
    * @fullName dexp(rate, logp)(x)
    * @memberof exponential
    */
   function dexp(rate, logp) {
      logp = logp === true;

      if (utils.hasNaN(rate) || rate < 0) { return function(x) { return NaN; }; }

      return function(x) {
         if (rate === Infinity) { return NaN; }
         if (x < 0) {
            return logp ? -Infinity : 0;
         }
         return logp ? -rate * x + Math.log(rate)
                     : rate * Math.exp(-rate * x);
      };
   }

   /**
    * Evaluates the lower-tail cdf at `x` for the exponential distribution:
    * $$\textrm{pexp}(\lambda)(x) = 1 - e^{-\lambda x}$$
    *
    * Expects the rate $\lambda > 0$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability instead.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    * @fullName pexp(rate, lowerTail, logp)(x)
    * @memberof exponential
    */
   function pexp(rate, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (utils.hasNaN(rate)) { return function(x) { return NaN; }; }

      return function(x) {
         var ret;

         if (utils.hasNaN(x)) { return NaN; }
         if (rate < 0 && rate !== -Infinity) { return NaN; }
         if (x <= 0) {
            ret = lowerTail ? 0 : 1;
            return logp ? Math.log(ret) : ret;
         }

         // basic behavior
         if (lowerTail) {
            return logp ? log1p(-Math.exp(-rate * x))
                        : -expm1(-rate * x);
         }
         // else right tail
         return logp ? -rate * x
                     : Math.exp(-rate * x);
      };
   }

   /**
    * Evaluates the exponential distribution's quantile function
    * (inverse cdf) at `p`:
    * $$\textrm{qexp}(\lambda)(p) = x \textrm{ such that } \textrm{prob}(X \leq x) = p$$
    * where $X$ is a random variable with the exponential distribution.
    *
    * Expects the rate $\lambda > 0$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
    * interpreted as an upper tail probability (returns
    * $x$ such that $\textrm{prob}(X > x) = p)$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
    * the logarithm of the desired probability.
    * @fullName qexp(rate, lowerTail, logp)(p)
    * @memberof exponential
    */
   function qexp(rate, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (utils.hasNaN(rate)) { return function(x) { return NaN; }; }

      return utils.qhelper(lowerTail, logp, 0, Infinity, function(p) {
         if (rate < 0 && rate !== -Infinity) { return NaN; }
         if (logp) {
            return lowerTail ? -Math.log(-expm1(p)) / rate
                             : -p / rate;
         }
         return lowerTail ? -log1p(-p) / rate
                          : -Math.log(p) / rate;
      });
   }

   // helper function for rexp
   sexp = (function() {
      var q, i, r;

      q = [Math.log(2)];
      r = 1;  // r will accumulate log(2)^k / k!
      for (i = 1; i < 16; i += 1) {
         r *= Math.log(2) / (i + 1);
         q[i] = q[i - 1] + r;
      }

      return function() {
         var a, u, ustar, umin;

         a = 0;
         do { u = rand(); } while (u <= 0 || u >= 1);
         u += u;
         while (u <= 1) {
            u += u;
            a += q[0];
         }
         u -= 1;

         if (u <= q[0]) { return a + u; }

         i = 0;
         umin = rand();
         do {
            ustar = rand();
            if (ustar < umin) { umin = ustar; }
            i += 1;
         } while (u > q[i]);

         return a + umin * q[0];
      };
   }());

   /**
    * Returns a random variate from the exponential distribution.
    *
    * Expects the rate $\lambda > 0$.
    *
    * Based on the R code in sexp.c (nmath)
    *
    * @fullName rexp(rate)()
    * @memberof exponential
    */
   function rexp(rate) {
      if (rate < 0) { return function() { return NaN; }; }

      return function() { return sexp() / rate; };
   }

   return {
      /**
       * Returns an object representing an exponential distribution with a
       * given rate $\lambda > 0$,
       * with properties `d`, `p`, `q`, `r`.
       * ```
       * expdistr(rate).d(x, logp)            // same as dexp(rate, logp)(x)
       * expdistr(rate).p(x, lowerTail, logp) // same as pexp(rate, lowerTail, logp)(x)
       * expdistr(rate).q(x, lowerTail, logp) // same as qexp(rate, lowerTail, logp)(x)
       * expdistr(rate).r()                   // same as rexp(rate)()
       * ```
       * @memberof exponential
       */
      expdistr: function(rate) {
         return {
            d: function(x, logp) { return dexp(rate, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pexp(rate, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qexp(rate, lowerTail, logp)(p);
            },
            r: function() { return rexp(rate)(); }
         };
      },
      dexp: dexp,
      pexp: pexp,
      qexp: qexp,
      rexp: rexp
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
