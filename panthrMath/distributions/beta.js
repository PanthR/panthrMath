(function(define) {'use strict';
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
   var bratio, lbeta, solve, stirlerr, bd0, log1p, C, rgen;

   C = require('../constants');
   bd0 = require('../basicFunc/bd0').bd0;
   stirlerr = require('../basicFunc/stirlerr').stirlerr;
   log1p = require('../basicFunc/log1p').log1p;
   lbeta = require('../basicFunc/lbeta').lbeta;
   bratio = require('../basicFunc/bratio').bratio;
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

      return function(x) {
         var lb, p, n;
         if (x < 0 || x > 1) {
            lb = -Infinity;
         } else if (a <= 2 || b <= 2) {
            lb = (a - 1) * Math.log(x) + (b - 1) * log1p(-x) - lbeta(a, b);
         } else {
            p = x;
            x = a - 1;
            n = a + b - 2;
            lb = Math.log(a + b - 1) +
               stirlerr(n) - stirlerr(x) - stirlerr(n - x) -
               bd0(x, n * p) - bd0(n - x, n * (1 - p)) +
               0.5 * Math.log(n / (C.twopi * x * (n - x)));
         }
         return logp ? lb : Math.exp(lb);
      };
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

      return function(x) {
         var lp;
         if (a <= 0 || b <= 0) { return NaN; }
         if (x > 0 && x < 1) {
            return bratio(a, b, x, lowerTail, logp);
         }
         lp = x <= 0 ? 0 : 1;
         lp = lowerTail ? lp : 1 - lp;
         return logp ? Math.log(lp) : lp;
      };
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

      /* eslint-disable complexity */
      return function(p) {
         // need reasonable endpoints
         var l, u, incr, realLogP;

         if (a <= 0 || b <= 0) { return NaN; }

         if (!logp && !(p >= 0 && p <= 1)) { return NaN; }
         realLogP = logp ? p : Math.log(p);
         realLogP = lowerTail ? realLogP : -realLogP;

         if (realLogP === -Infinity) { return 0; }
         if (realLogP === +Infinity) { return 1; }

         incr = 1e-1;
         l = 1e-1;
         while (f(l) > p === lowerTail) {
            l = l * incr;
            if (l < 1e-80) {
               return new Error('qbeta failed to find left endpoint');
            }
         }
         u = 1 - incr;
         while (f(u) < p === lowerTail) {
            incr *= .1;
            u = 1 - incr;
            if (incr < 1e-80) {
               return new Error('qbeta failed to find right endpoint');
            }
         }
         /* eslint-enable complexity */

         return solve(function(x) { return f(x); }, p, l, u);
      };
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
         };
         /* eslint-enable complexity */
      }

      // Algorithm BB
      beta = Math.sqrt((alpha - 2) / (2 * a0 * b0 - alpha));
      gamma = a0 + 1 / beta;

      return function() {
         var u1, u2, v, w, z, r, s;

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
         return a === a0 ? b0 / (b0 + w) : w / (b0 + w);
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
