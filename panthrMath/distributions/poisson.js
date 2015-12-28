(function(define) {'use strict';
define(function(require) {

   // No input validation provided.

   /**
    * Provides probability mass function, cumulative distribution function,
    * quantile function, and random number generator
    * for the Poisson distribution, which is defined by the pmf
    * $$f(x;\lambda) = \frac{\lambda^x}{x!}e^{-\lambda}$$
    * where $x=0,1,2,\ldots$ and $\lambda > 0$ is the mean.
    *
    * `dpois` provides access to this probability mass function,
    * `ppois` to the cumulative distribution function, `qpois` to the
    * quantile function (inverse cdf)
    * and `rpois` to random deviates.
    *
    * Finally, you can use `pois` to obtain an object
    * representing the Poisson distribution for a given value of $\lambda$.
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

   /**
    * Evaluates the Poisson pmf at `x`:
    * $$\textrm{dpois}(\lambda)(x) = \frac{\lambda^x}{x!}e^{-\lambda}$$
    *
    * Expects $\lambda > 0$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the
    * logarithm of the result.
    * @fullName dpois(lambda, logp)(x)
    * @memberof poisson
    */
   function dpois(lambda, logp) {
      logp = logp === true;
      return function(x) {
         return logp ? lpoisson(lambda)(x) : Math.exp(lpoisson(lambda)(x));
      };
   }

   /**
    * Evaluates the lower-tail cdf at `x` for the Poisson distribution:
    * $$\textrm{ppois}(\lambda)(x) = e^{-\lambda} \sum_{i=0}^{\left \lfloor{x}\right \rfloor}\frac{\lambda^i}{i!}$$
    *
    * Expects $\lambda > 0$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability ($P(X > x)$) instead.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    *
    * @fullName ppois(lambda, lowerTail, logp)(x)
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
    * Evaluates the Poisson distribution's quantile function
    * (inverse cdf) at `p`:
    * $$\textrm{qpois}(\lambda)(p) = x \textrm{ such that } \textrm{prob}(X \leq x) = p$$
    * where $X$ is a random variable with the $\textrm{Poisson}(\lambda)$ distribution.
    *
    * Expects $\lambda > 0$ and $0 \leq p \leq 1$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
    * interpreted as an upper tail probability (returns
    * $x$ such that $\textrm{prob}(X > x) = p)$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
    * the logarithm of the desired probability.
    *
    * @fullName qpois(lambda, lowerTail, logp)(p)
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
    * Returns a random variate from the $\textrm{Poisson}(\lambda)$ distribution.
    *
    * Expects $\lambda > 0$.
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
       * Returns an object representing a Poisson distribution
       * for a given value of $\lambda > 0$, with properties `d`, `p`, `q`, `r`.
       * ```
       * pois(a, b).d(x, logp)            // same as dpois(a, b, logp)(x)
       * pois(a, b).p(x, lowerTail, logp) // same as ppois(a, b, lowerTail, logp)(x)
       * pois(a, b).q(x, lowerTail, logp) // same as qpois(a, b, lowerTail, logp)(x)
       * pois(a, b).r()                   // same as rpois(a, b)()
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
            r: function() { return rpois(lambda)(); }
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
