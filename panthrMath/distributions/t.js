(function(define) {'use strict';
define(function(require) {

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the Student's t distribution, which is defined by the pdf
    * $$p(x;\nu) = \frac{1}{\sqrt{\nu}\,\textrm{B}(1/2, \nu/2)} \left(1+\frac{x^2}{\nu} \right)^{-\frac{\nu + 1}{2}}$$
    * where B is the Beta function, $x \in (-\infty, \infty)$ and $\nu >0$ is the *degrees of freedom*.
    *
    * `dt` provides access to this probability density function,
    * `pt` to the cumulative distribution function, `qt` to the
    * quantile function (inverse cdf)
    * and `rt` to random deviates.
    *
    * Finally, you can use `tdistr` to obtain an object
    * representing the t distribution for a given value of $\nu$.
    *
    * @module distributions.t
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var C, bd0, stirlerr, sqrt2pi, bratio, log1p, solve, rgen;
   // var expm1,  qnorm;

   C = require('../constants');
   bd0 = require('../basicFunc/bd0').bd0;
   stirlerr = require('../basicFunc/stirlerr').stirlerr;
   sqrt2pi = C.sqrt2pi;
   bratio = require('../basicFunc/bratio').bratio;
   log1p = require('../basicFunc/log1p').log1p;
   // expm1 = require('../basicFunc/expm1').expm1;
   solve = require('../utils').binSearchSolve;
   // qnorm = require('./normal').qnorm;
   rgen = require('../rgen/rgen');

   /*
    * Return the log-of-density function for student's t distribution.
    * `df` is the degrees of freedom.
    * See: R source code.
    */
   function dtlog(df) {
      var t;
      t = -bd0(df / 2, (df + 1) / 2) + stirlerr((df + 1) / 2) - stirlerr(df / 2);
      return function(x) {
         var x2n, logx2n, u;
         x2n = x * x / df;
         logx2n = Math.log(1 + x2n) / 2;
         u = x2n > 0.2 ? df * logx2n
                       : x * x / 2 - bd0(df / 2, (df + x * x) / 2);
         return t - u - (logx2n + Math.log(sqrt2pi));
      };
   }

   // density function
   // dt(df) returns a function for calculating cumulative distribution
   /**
    * Evaluates the t density function at `x`:
    * $$\textrm{dt}(\textrm{df})(x) = \frac{1}{\sqrt{\textrm{df}}\,\textrm{B}(1/2, \textrm{df}/2)} \left(1+\frac{x^2}{\textrm{df}} \right)^{-\frac{\textrm{df} + 1}{2}}$$
    *
    * Expects $\textrm{df} > 0$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the
    * logarithm of the result.
    *
    * @fullName dpois(lambda, logp)(x)
    * @fullName dt(df, logp)(x)
    * @memberof t
    */
   function dt(df, logp) {
      logp = logp === true;
      var dtl;
      dtl = dtlog(df);
      return function(x) {
         return logp ? dtl(x) : Math.exp(dtl(x));
      };
   }

   // cumulative distribution function, log version
   // based on pt.c from R implementation
   function ptlog(df, lowerTail) {
      return function(x) {
         var val;

         val = df > x * x ?
               bratio(0.5, df / 2, x * x / (df + x * x), false, true) :
               bratio(df / 2, 0.5, 1 / (1 + x / df * x), true, true);

         if (x <= 0) { lowerTail = !lowerTail; }

         return lowerTail ? log1p(-0.5 * Math.exp(val))
                          : val - Math.log(2);
      };
   }

   // cumulative distribution function
   /**
    * Evaluates the lower-tail cdf at `x` for the t distribution:
    * $$\textrm{pt}(\textrm{df})(x) = \int_{-\infty}^{x} \frac{1}{\sqrt{\textrm{df}}\,\textrm{B}(1/2, \textrm{df}/2)} \left(1+\frac{x^2}{\textrm{df}} \right)^{-\frac{\textrm{df} + 1}{2}}$$
    *
    * Expects $\textrm{df} > 0$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability instead.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    * @fullName pt(df, lowerTail, logp)(x)
    * @memberof t
    */
   function pt(df, lowerTail, logp) {
      var ptl;

      lowerTail = lowerTail !== false;
      logp = logp === true;
      ptl = ptlog(df, lowerTail);
      return function(x) {
         return logp ? ptl(x) : Math.exp(ptl(x));
      };
   }

   // inverse cdf
   // From qt.c in R
   /**
    * Evaluates the t distribution's quantile function
    * (inverse cdf) at `p`:
    * $$\textrm{qt}(\textrm{df})(p) = x \textrm{ such that } \textrm{prob}(X \leq x) = p$$
    * where $X$ is a random variable with the $t(\textrm{df})$ distribution.
    *
    * Expects $\textrm{df} > 0$ and $0 \leq p \leq 1$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
    * interpreted as an upper tail probability (returns
    * $x$ such that $\textrm{prob}(X > x) = p)$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
    * the logarithm of the desired probability.
    *
    * Based on: R code
    *
    * @fullName qt(df, lowerTail, logp)(p)
    * @memberof t
    */
   function qt(df, lowerTail, logp) {
      lowerTail = lowerTail !== false;
      logp = logp === true;

      return function(p) {
         // var pp;

         // two-tailed prob, pp = 2 * min(p, 1-p)
         // pp = 2 * (logp ? Math.min(Math.exp(p), -expm1(p))
                        // : Math.min(p, 1 - p));

         // Using the solver on whole range (possibly inefficient but works)
         // if (df < 1) {
            return solve(function(x) {
               return pt(df, lowerTail, logp)(x);
            }, p);
         // }

         // if (df > 1e20) {
         //    return qnorm(0, 1, lowerTail, logp)(p);
         // }

         // TODO: Could consider optimizing for df close to 1 or 2

      };
   }

   /**
    * Returns a random variate from the $t(\textrm{df})$ distribution.
    *
    * Expects $\textrm{df} > 0$.
    * @memberof t
    */
   function rt(df) {
         return function() {
            var u1, u2, rsq;
            do {
               u1 = rgen.random() * 2 - 1;
               u2 = rgen.random() * 2 - 1;
               rsq = u1 * u1 + u2 * u2;
            } while (rsq >= 1);
            return u1 * Math.sqrt(df * ( Math.pow(rsq, -2 / df) - 1) / rsq);
         };
      }

   return {
      /**
       * Returns an object representing a t distribution with $\textrm{df} > 0$
       * degrees of freedom, with properties `d`, `p`, `q`, `r`.
       * ```
       * tdistr(a, b).d(x, logp)            // same as dt(a, b, logp)(x)
       * tdistr(a, b).p(x, lowerTail, logp) // same as pt(a, b, lowerTail, logp)(x)
       * tdistr(a, b).q(x, lowerTail, logp) // same as qt(a, b, lowerTail, logp)(x)
       * tdistr(a, b).r()                   // same as rt(a, b)()
       * ```
       * @memberof t
       */
      tdistr: function(df) {
         return {
            d: function(x, logp) { return dt(df, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pt(df, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qt(df, lowerTail, logp)(p);
            },
            r: function() { return rt(df)(); }
         };
      },
      dt: dt,
      pt: pt,
      qt: qt,
      rt: rt
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
