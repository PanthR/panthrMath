(function(define) {
'use strict';
define(function(require) {

   /**
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
   var rgen, utils;

   rgen = require('../rgen/rgen');
   utils = require('../utils');

   /**
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
   function dunif(a, b, logp) {
      logp = logp === true;

      if (utils.hasNaN(a, b) || b <= a) {
         return function(x) { return NaN; };
      }

      return function(x) {
         var p;

         if (utils.hasNaN(x)) { return NaN; }
         p = x < a || x > b ? 0 : 1 / (b - a);

         return logp ? Math.log(p) : p;
      };
   }

   /**
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
   function punif(a, b, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (b < a || !utils.isFinite(a) || !utils.isFinite(b)) {
         return function(x) { return NaN; };
      }

      return function(x) {
         if (utils.hasNaN(x)) { return NaN; }
         if (x >= b) { return utils.adjustLower(1, lowerTail, logp); }
         if (x <= a) { return utils.adjustLower(0, lowerTail, logp); }

         return utils.adjustLower(
            (lowerTail ? x - a : b - x) / (b - a), true, logp
         );
      };
   }

   /**
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
   function qunif(a, b, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (b < a || !utils.isFinite(a) || !utils.isFinite(b)) {
         return function(x) { return NaN; };
      }

      return function(p) {
         p = logp ? Math.exp(p) : p;

         if (p < 0 || p > 1) { return NaN; }
         return lowerTail ? a + p * (b - a) : b - p * (b - a);
      };
   }

   /**
    * Returns a random variate from the $\textrm{Uniform}(a,b)$ distribution.
    *
    * Expects $a<b$.
    * @fullName runif(a, b)()
    * @memberof uniform
    */
   function runif(a, b) {
      return function() {
         return a + (b - a) * rgen.random();
      };
   }

   return {
      /**
       * Returns an object representing a Uniform distribution on $[a, b]$,
       * with properties `d`, `p`, `q`, `r`.
       * ```
       * unif(a, b).d(x, logp)            // same as dunif(a, b, logp)(x)
       * unif(a, b).p(x, lowerTail, logp) // same as punif(a, b, lowerTail, logp)(x)
       * unif(a, b).q(x, lowerTail, logp) // same as qunif(a, b, lowerTail, logp)(x)
       * unif(a, b).r()                   // same as runif(a, b)()
       * ```
       * @memberof uniform
       */
      unif: function(a, b) {
         return {
            d: function(x, logp) { return dunif(a, b, logp)(x); },
            p: function(q, lowerTail, logp) {
               return punif(a, b, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qunif(a, b, lowerTail, logp)(p);
            },
            r: function() { return runif(a, b)(); }
         };
      },
      dunif: dunif,
      punif: punif,
      qunif: qunif,
      runif: runif
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
