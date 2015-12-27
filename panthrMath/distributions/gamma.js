(function(define) {'use strict';
define(function(require) {

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the Gamma distribution, which is defined by the pdf
    * $$p(x;a,s) = \frac{1}{s^a\Gamma(a)}x^{a-1}e^{-x/s}$$
    * where $x\in(0,\infty)$. The parameter $a>0$ is called the *shape*
    * and the parameter $s>0$ the *scale*.
    *
    * `dgamma` provides access to this probability density function,
    * `pgamma` to the cumulative distribution function, `qgamma` to the
    * quantile function (inverse cdf)
    * and `rgamma` to random deviates.
    *
    * Finally, you can use `gammadistr` to obtain an object
    * representing the distribution for some values of the parameters.
    * @module distributions.gamma
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var lpoisson, gratio, exponential, rgen;

   lpoisson = require('../basicFunc/lpoisson').lpoisson;
   gratio = require('../basicFunc/gratio');
   exponential = require('../rgen/exponential');
   rgen = require('../rgen/rgen');


   // helper function
   function dgammaLog(a, s) {
      if (a < 1) {
         return function(x) {
            if (x === 0) { return Infinity; }
            return Math.log(a / x) + lpoisson(x / s)(a);
         };
      }
      if (a === 1) {
         return function(x) {
            return -x / s - Math.log(s);
         };
      }
      // a > 1
      return function(x) {
         if (x === 0) { return -Infinity; }
         return -Math.log(s) + lpoisson(x / s)(a - 1);
      };
   }

   // density function
   // 1 / (s^a * gamma(a)) * x ^ (a - 1) * exp(-x / s)
   /**
    * Evaluates the Gamma density function at `x`:
    * $$\textrm{dgamma}(a, s)(x) = \frac{1}{s^a\Gamma(a)}x^{a-1}e^{-x/s}$$
    *
    * Expects $a > 0$, $s > 0$, and $x > 0$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the
    * logarithm of the result.
    *
    * @fullName dgamma(a, s, logp)(x)
    * @memberof gamma
    */
   function dgamma(a, s, logp) {
      logp = logp === true;
      return function(x) {
         return logp ? dgammaLog(a, s)(x) : Math.exp(dgammaLog(a, s)(x));
      };
   }

   // cumulative distribution function
   // a = shape, s = scale (rate = 1 / s)
   /**
    * Evaluates the Gamma cumulative distribution
    * function at `x` (lower tail probability):
    * $$\textrm{pgamma}(a, s)(x) = \frac{1}{s^a\Gamma(a)}\int_0^x t^{a-1}e^{-t/s}dt$$
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability instead.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    *
    * Expects $a > 0$, $s > 0$, and $x > 0$.
    *
    * Based on: *Computation of the Incomplete Gamma Function Ratios and their
    * Inverse*, by DiDonato and Morris, 1992
    * @fullName pgamma(a, s, lowerTail, logp)(x)
    * @memberof gamma
    */
   function pgamma(a, s, lowerTail, logp) {
      var f;

      if (s <= 0 || a <= 0) { return function(x) { return NaN; }; }
      logp = logp === true;
      lowerTail = lowerTail !== false;
      f = lowerTail ? gratio.gratio(a) : gratio.gratioc(a);

      return function(x) {
         var p;

         p = x > 0 ? f(x / s) : lowerTail ? 0 : 1;
         return logp ? Math.log(p) : p;
      };
   }

   // inverse cdf
   // a = shape, s = scale (rate = 1 / s)
   /**
    * Evaluates the Gamma distribution's quantile function (inverse cdf) at `p`:
    * $$\textrm{qgamma}(a, s)(p) = x \textrm{ such that } \textrm{prob}(X \leq x) = p$$
    * where $X$ is a random variable with the $\textrm{Gamma}(a, s)$ distribution.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
    * interpreted as an upper tail probability (returns
    * $x$ such that $\textrm{prob}(X > x) = p)$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
    * the logarithm of the desired probability.
    *
    * Expects $a>0$, $s>0$, and $0 \leq p \leq 1$.
    * @fullName qgamma(a, s, lowerTail, logp)(p)
    * @memberof gamma
    */
   function qgamma(a, s, lowerTail, logp) {
      var gaminv;

      if (s <= 0 || a <= 0) { return function(x) { return NaN; }; }
      logp = logp === true;
      lowerTail = lowerTail !== false;
      gaminv = gratio.gaminv(a); // function of p, lowerTail

      return function(p) {
         p = logp ? Math.exp(p) : p;
         if (!(p >= 0 && p <= 1)) { return NaN; }
         return s * gaminv(p, lowerTail);
      };
   }

   // rgamma stuff...
   // Client calls one of three methods depending on the value of a.
   // a for 'shape', s for 'scale'
   /**
    * Returns a random variate from the $\textrm{Gamma}(a, s)$ distribution.
    *
    * Expects $a > 0$ and $s > 0$.
    *
    * TODO: Add links to algorithms
    * @memberof gamma
    */
   function rgamma(a, s) {
      if (a < 1) {
         return gammaBest(a, s);
      }
      if (a > 1) {
         return gammaCheng(a, s);
      }
      return exponential(1 / s);
   }

   // a < 1:  Best/Ahrens/Dieter Algorithm
   function gammaBest(a, s) {
      var t, b;
      t = .07 + .75 * Math.sqrt(1 - a);
      b = 1 + Math.exp(-t) * a / t;

      return function() {
         var x, y, u1, u2, v;
         /* eslint-disable no-constant-condition */
         while (true) {
         /* eslint-enable */
            u1 = rgen.random();
            u2 = rgen.random();
            v = b * u1;
            if (v <= 1) {
               x = t * Math.pow(v, 1 / a);
               if (u2 <= (2 - x) / (2 + x) || u2 <= Math.exp(-x)) {
                  return x * s;
               }
            } else {
               x = -Math.log(t * (b - v) / a);
               y = x / t;
               if (u2 * (a + y * (1 - a)) <= 1 ||
                   u2 <= Math.pow(y, a - 1))
               {
                  return x * s;
               }
            }
         }
      };
   }

   // a > 1:  Cheng/Feast Algorithm
   function gammaCheng(a, s) {
      var c1, c2, c3, c4;
      c1 = a - 1;
      c2 = (a - 1 / (6 * a)) / c1;
      c3 = 2 / c1;
      c4 = c3 + 2;

      return function() {
         var u1, u2, w;
         /* eslint-disable no-constant-condition */
         while (true) {
         /* eslint-enable */
            u1 = rgen.random();
            u2 = rgen.random();
            w = c2 * u2 / u1;
            if (c3 * u1 + w + 1 / w <= c4) {
               return c1 * w * s;
            }
            if (c3 * Math.log(u1) - Math.log(w) + w < 1) {
               return c1 * w * s;
            }
         }
      };
   }

   return {
      /**
       * Returns an object representing a gamma distribution, with properties `d`, `p`, `q`, `r`.
       * ```
       * gammadistr(a, s).d(x, logp)            // same as dgamma(a, s, logp)(x)
       * gammadistr(a, s).p(x, lowerTail, logp) // same as pgamma(a, s, lowerTail, logp)(x)
       * gammadistr(a, s).q(x, lowerTail, logp) // same as qgamma(a, s, lowerTail, logp)(x)
       * gammadistr(a, s).r()                   // same as rgamma(a, s)()
       * ```
       * @memberof gamma
       */
      gammadistr: function(a, s) {
         return {
            d: function(x, logp) { return dgamma(a, s, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pgamma(a, s, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qgamma(a, s, lowerTail, logp)(p);
            },
            r: function() { return rgamma(a, s)(); }
         };
      },
      dgamma: dgamma,
      pgamma: pgamma,
      qgamma: qgamma,
      rgamma: rgamma
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
