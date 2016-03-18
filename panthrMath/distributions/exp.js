(function(define) {
'use strict';
define(function(require) {

   /** TODO
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the continuous uniform distribution on $[a, b]$, which
    * is defined by the pdf
    * $$ f(x; a, b) = \begin{cases}
    *   \frac{1}{b-a},  & \text{if $a \leq x \leq b$} \\\\
    *   0, & \text{if $x < a$ or $x > b$} \end{cases} $$
    * where $a < b$ and $x \in (-\infty, \infty)$.
    *
    * `dunif` provides access to this probability density function,
    * `punif` to the cumulative distribution function, `qunif` to the
    * quantile function (inverse cdf)
    * and `runif` to random deviates.
    *
    * Finally, you can use `unif` to obtain an object
    * representing the Uniform distribution for given values of
    * `a` and `b`.
    *
    * @module distributions.uniform
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var sexp, log1p, expm1, rand;

   log1p = require('../basicFunc/log1p').log1p;
   expm1 = require('../basicFunc/expm1').expm1;
   rand = require('../rgen/rgen').random;

   /** TODO
    * Evaluates the Uniform density function at `x`:
    * $$\textrm{dunif}(a, b)(x) = \begin{cases}
    *   \frac{1}{b-a},  & \text{if $a \leq x \leq b$} \\\\
    *   0, & \text{if $x < a$ or $x > b$} \end{cases} $$
    *
    * Expects $a < b$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the
    * logarithm of the result.
    *
    * @fullName dunif(a, b, logp)(x)
    * @memberof uniform
    */
   function dexp(rate, logp) {
      logp = logp === true;

      if (rate <= 0) { return function(x) { return NaN; }; }

      return function(x) {
         if (x < 0) {
            return logp ? -Infinity : 0;
         }
         return logp ? -rate * x + Math.log(rate)
                     : rate * Math.exp(-rate * x);
      };
   }

   /** TODO
    * Evaluates the lower-tail cdf at `x` for the Uniform distribution:
    * $$\textrm{punif}(a, b)(x) = \begin{cases}
    *   \frac{x-a}{b-a},  & \text{if $a \leq x \leq b$} \\\\
    *   0,                & \text{if $x < a$} \\\\
    *   1,                & \text{if $x > b$} \end{cases} $$
    *
    * Expects $a < b$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability instead.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    * @fullName punif(a, b, lowerTail, logp)(x)
    * @memberof uniform
    */
   function pexp(rate, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (rate <= 0) { return function(x) { return NaN; }; }

      return function(x) {
         var ret;

         if (x < 0) {
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

   /** TODO
    * Evaluates the Uniform distribution's quantile function
    * (inverse cdf) at `p`:
    * $$\textrm{qunif}(a, b)(p) = x \textrm{ such that } \textrm{prob}(X \leq x) = p$$
    * where $X$ is a random variable with the $\textrm{Uniform}(a, b)$ distribution.
    *
    * Expects $a < b$ and $0 \leq p \leq 1$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
    * interpreted as an upper tail probability (returns
    * $x$ such that $\textrm{prob}(X > x) = p)$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
    * the logarithm of the desired probability.
    * @fullName qunif(a, b, lowerTail, logp)(p)
    * @memberof uniform
    */
   function qexp(rate, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (rate <= 0) { return function(x) { return NaN; }; }

      return function(p) {
         if (logp) {
            if (p > 0) { return NaN; }
            return lowerTail ? -Math.log(-expm1(p)) / rate
                             : -p / rate;
         }
         // else, not log
         if (!(p >= 0 && p <= 1)) { return NaN; }
         return lowerTail ? -log1p(-p) / rate
                          : -Math.log(p) / rate;
      };
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

   /** TODO
    * See: R code, sexp.c (nmath)
    *
    *
    * @memberof exp
    */
   function rexp(rate) {
      if (rate < 0) { return function() { return NaN; }; }

      return function() { return rate * sexp(); };
   }

   return {
      /** TODO
       * Returns an object representing an exponential distribution with a
       * given rate (`rate` > 0),
       * with properties `d`, `p`, `q`, `r`.
       * ```
       * expdistr(rate).d(x, logp)            // same as dexp(rate, logp)(x)
       * expdistr(rate).p(x, lowerTail, logp) // same as pexp(rate, lowerTail, logp)(x)
       * expdistr(rate).q(x, lowerTail, logp) // same as qexp(rate, lowerTail, logp)(x)
       * expdistr(rate).r()                   // same as rexp(rate)()
       * ```
       * @memberof exp
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
