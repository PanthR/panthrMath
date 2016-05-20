(function(define) {
'use strict';
define(function(require) {

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the Beta distribution, which is defined by the pdf
    * $$f(x;a,b) = \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}x^{(a-1)}(1-x)^{(b-1)}$$
    * where $x\in[0,1]$ and the parameters $a,b > 0$.
    *
    * `dbeta` provides access to this probability density function,
    * `pbeta` to the cumulative distribution function, `qbeta` to the
    * quantile function (inverse cdf)
    * and `rbeta` to random deviates.
    *
    * Finally, you can use `betadistr` to obtain an object
    * representing the distribution for some values of the parameters.
    * @module distributions.beta
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var bratio, lbeta, solve, log1p, rgen, dbinomLog, utils;

   utils = require('../utils');
   log1p = require('../basicFunc/log1p').log1p;
   lbeta = require('../basicFunc/lbeta').lbeta;
   bratio = require('../basicFunc/bratio').bratio;
   dbinomLog = require('../basicFunc/dbinomlog').dbinomLog;
   solve = require('../utils').binSearchSolve;
   rgen = require('../rgen/rgen');

   /**
    * Evaluates the Beta density function at `x`:
    * $$\textrm{dbeta}(a,b)(x) = \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}x^{(a-1)}(1-x)^{(b-1)}$$
    *
    * Expects $a > 0$, $b > 0$, and $0 \leq x \leq 1$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the
    * logarithm of the result.
    *
    * @fullName dbeta(a, b, logp)(x)
    * @memberof beta
    */
   function dbeta(a, b, logp) {
      logp = logp === true;

      /* eslint-disable complexity */
      return function(x) {
         var lb, p, n, alogx, blogmx;

         if (utils.hasNaN(x, a, b) || a < 0 || b < 0) { return NaN; }
         if (x < 0 || x > 1) {
            lb = -Infinity;
         } else if (a === 0 && b === 0) { // point mass 1/2 at each of {0,1} :
            lb = x === 0 || x === 1 ? Infinity : -Infinity;
         } else if (a === Infinity && b === Infinity) {
            lb = x === 0.5 ? Infinity : -Infinity; // point mass 1 at 1/2
         } else if (a / b === 0) { // a=0 or b=Inf. point mass 1 at 0
            lb = x === 0 ? Infinity : -Infinity;
         } else if (b / a === 0) { // b=0 or a=Inf. point mass 1 at 1
            lb = x === 1 ? Infinity : -Infinity;
         } else if (a <= 2 || b <= 2) {
            // Direct beta calculation
            alogx = a === 1 && x === 0 ? 0 : (a - 1) * Math.log(x);
            blogmx = b === 1 && x === 1 ? 0 : (b - 1) * log1p(-x);
            lb = alogx + blogmx - lbeta(a, b);
         } else {
            p = x;
            x = a - 1;
            n = a + b - 2;
            lb = Math.log(a + b - 1) + dbinomLog(n, p)(x);
         }
         return logp ? lb : Math.exp(lb);
      };
      /* eslint-enable complexity */
   }

   /**
    * Evaluates the Beta cumulative distribution
    * function at `x` (lower tail probability):
    * $$\textrm{pbeta}(a, b)(x) = I_x(a, b)=G(a,b)\int_0^xt^{a-1}(1-t)^{b-1}dt$$
    * where $B(a,b)=1/G(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b)$ is the
    * usual beta function.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability instead.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    *
    * Expects $a>0$, $b>0$, and $0 \leq x \leq 1$.
    *
    * Based on: *Algorithm 708: Significant Digit Computation of the Incomplete Beta Function
    * Ratios*, by DiDonato and Morris, 1992
    * @fullName pbeta(a, b, lowerTail, logp)(x)
    * @memberof beta
    */
   function pbeta(a, b, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      /* eslint-disable complexity */
      return function(x) {
         if (utils.hasNaN(x, a, b) || a < 0 || b < 0) { return NaN; }

         if (x <= 0) { return utils.adjustLower(0, lowerTail, logp); }
         if (x >= 1) { return utils.adjustLower(1, lowerTail, logp); }

         // Point mass limit cases
         if (a === 0 && b === 0) { // point mass 1/2 at each of {0,1}
            return logp ? Math.log(0.5) : 0.5;
         }
         if (a === Infinity && b === Infinity) { // point mass 1 at 1/2
            return utils.adjustLower(x < 0.5 ? 0 : 1, lowerTail, logp);
         }
         if (a / b === 0) { // point mass 1 at 0 ==> P(X <= x) = 1, all x > 0
            return utils.adjustLower(1, lowerTail, logp);
         }
         if (b / a === 0) { // point mass 1 at 1 ==> P(X <= x) = 0, all x < 1
            return utils.adjustLower(0, lowerTail, logp);
         }
         return bratio(a, b, x, lowerTail, logp);
      };
      /* eslint-enable complexity */
   }

   // inverse cdf
   // preliminary implementation uses binSearchSolve, similar to our
   /**
    * Evaluates the Beta quantile function (inverse cdf) at `p`:
    * $$\textrm{qbeta}(a, b)(p) = x \textrm{ such that } \textrm{prob}(X \leq x) = p$$
    * where $X$ is a random variable with the $\textrm{Beta}(a,b)$ distribution.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
    * interpreted as an upper tail probability (returns
    * $x$ such that $\textrm{prob}(X > x) = p)$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
    * the logarithm of the desired probability.
    *
    * Expects $a>0$, $b>0$, and $0 \leq p \leq 1$.
    * @fullName qbeta(a, b, lowerTail, logp)(p)
    * @memberof beta
    */
   function qbeta(a, b, lowerTail, logp) {
      var f;

      lowerTail = lowerTail !== false;
      logp = logp === true;
      f = pbeta(a, b, lowerTail, logp);

      if (utils.hasNaN(a, b) || a < 0 || b < 0) {
         return function() { return NaN; };
      }

      /* eslint-disable complexity */
      return utils.qhelper(lowerTail, logp, 0, 1, function(p) {
         // need reasonable endpoints
         var l, u, incr, unloggedp;

         if (utils.hasNaN(p)) { return NaN; }

         unloggedp = logp ? Math.exp(p) : p;

         if (a === 0 && b === 0) { // point mass 1/2 at each of {0,1}
            // WARNING: Following R here, and ignoring the question of
            // upper/lower tail.
            return unloggedp < 0.5 ? 0
                 : unloggedp > 0.5 ? 1
                 : 0.5;
         }
         if (a === Infinity && b === Infinity) { // point mass 1 at 1/2
            return 0.5;
         }
         if (a / b === 0) { // point mass 1 at 0 ==> P(X <= x) = 1, all x > 0
            return 0;
         }
         if (b / a === 0) { // point mass 1 at 1 ==> P(X <= x) = 0, all x < 1
            return 1;
         }

         incr = 1e-1;
         l = 1e-1;
         while (f(l) > p === lowerTail) {
            l *= incr;
            if (l < 1e-80) {
               return new Error('qbeta failed to find left endpoint');
            }
         }
         u = 1 - incr;
         while (f(u) < p === lowerTail) {
            incr *= 0.1;
            u = 1 - incr;
            if (incr < 1e-80) {
               return new Error('qbeta failed to find right endpoint');
            }
         }
         /* eslint-enable complexity */

         return solve(function(x) { return f(x); }, p, l, u);
      });
   }

   // Following R's code, which follows:
   // Cheng 1978 "Generating beta variates with nonintegral shape parameters"
   /**
    * Returns a random variate from the $\textrm{Beta}(a, b)$ distribution.
    *
    * Expects $a>0$ and $b>0$.
    *
    * Based on R's code; see: *Generating beta variates with nonintegral shape parameters*, by
    * RCH Cheng, 1978
    * @fullName rbeta(a, b)()
    * @memberof beta
    */
   function rbeta(a, b) {
      var a0, b0, alpha, beta, gamma, delta, k1, k2;

      if (a <= 0 || b <= 0) { return function() { return NaN; }; }

      a0 = Math.min(a, b);
      b0 = Math.max(a, b);
      alpha = a0 + b0;

      if (a0 <= 1) {
         // Algorithm BC
         beta = 1 / a0;
         delta = 1 + b0 - a0;
         k1 = delta * (0.0138889 + 0.0416667 * a0) / (b0 * beta - 0.777778);
         k2 = 0.25 + (0.5 + 0.25 / delta) * a0;

         /* eslint-disable complexity */
         return function() {
            var u1, u2, v, w, y, z;

         /* eslint-disable no-constant-condition */
            while (true) {
               u1 = rgen.random();
               u2 = rgen.random();
               y = u1 * u2;
               z = u1 * y;
               v = beta * Math.log(u1 / (1 - u1));
               w = b0 * Math.exp(v);
               if (w === Infinity) { w = Number.MAX_VALUE; }
               if (u1 < 0.5) {
                  if (0.25 * u2 + z - y < k1 &&
                     alpha * (Math.log(alpha / (a0 + w)) + v) - 1.3862944 >= Math.log(z)) {
                     break;
                  }
               } else if (z <= 0.25 ||
                  z < k2 && alpha * (Math.log(alpha / (a0 + w)) + v) - 1.3862944 >= Math.log(z)) {
                  break;
               }
            }
            return a === a0 ? a0 / (a0 + w) : w / (a0 + w);
         /* eslint-enable no-constant-condition */
         };
         /* eslint-enable complexity */
      }

      // Algorithm BB
      beta = Math.sqrt((alpha - 2) / (2 * a0 * b0 - alpha));
      gamma = a0 + 1 / beta;

      return function() {
         var u1, u2, v, w, z, r, s;

         /* eslint-disable no-constant-condition */
         while (true) {
            u1 = rgen.random();
            u2 = rgen.random();

            v = beta * Math.log(u1 / (1 - u1));
            w = a0 * Math.exp(v);
            if (w === Infinity) { w = Number.MAX_VALUE; }

            z = u1 * u1 * u2;
            r = gamma * v - 1.3862944;
            s = a + r - w;
            if (s + 2.609438 >= 5 * z || s > Math.log(z) ||
               r + alpha * Math.log(alpha / (b0 + w)) >= Math.log(z)) {
               break;
            }
         }
         /* eslint-enable no-constant-condition */
         return a !== a0 ? b0 / (b0 + w) : w / (b0 + w);
      };
   }

   return {
      /**
       * Returns an object representing a beta distribution, with properties `d`, `p`, `q`, `r`.
       * ```
       * betadistr(a, b).d(x, logp)            // same as dbeta(a, b, logp)(x)
       * betadistr(a, b).p(x, lowerTail, logp) // same as pbeta(a, b, lowerTail, logp)(x)
       * betadistr(a, b).q(x, lowerTail, logp) // same as qbeta(a, b, lowerTail, logp)(x)
       * betadistr(a, b).r(n)                  // same as rbeta(a, b)(n)
       * ```
       * @memberof beta
       */
      betadistr: function(a, b) {
         return {
            d: function(x, logp) { return dbeta(a, b, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pbeta(a, b, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qbeta(a, b, lowerTail, logp)(p);
            },
            r: function(n) { return rbeta(a, b)(n); }
         };
      },
      dbeta: dbeta,
      pbeta: pbeta,
      qbeta: qbeta,
      rbeta: rbeta
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
