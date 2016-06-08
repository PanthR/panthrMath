(function(define) {
'use strict';
define(function(require) {

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the negative binomial distribution with parameters `size` and `prob`.
    *
    * This represents the number of failures which occur in a sequence of Bernoulli
    * trials before a target number (`size`) of successes is reached. The parameter
    * `prob` represents the probability of success in each trial.
    *
    * `size` must be strictly positive, and `prob` must be in `(0, 1]`.
    *
    * `dnbinom` provides access to the probability density function,
    * `pnbinom` to the cumulative distribution function, `qnbinom` to the
    * quantile function (inverse cdf)
    * and `rnbinom` to random deviates.
    *
    * Finally, you can use `nbinom` to obtain an object
    * representing the negative binomial distribution for given values of
    * `size` and `prob`.
    *
    * @module distributions.nbinom
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var beta, dbinomLog, adjustLower, qnorm, pWrap, discInvCdf, discInvCdf2, rpois, rgamma;

   dbinomLog = require('../basicFunc/dbinomLog').dbinomLog;
   beta = require('./beta');
   qnorm = require('./normal').qnorm;
   adjustLower = require('../utils').adjustLower;
   pWrap = require('../utils').pWrap;
   discInvCdf = require('../utils').discInvCdf;
   discInvCdf2 = require('../utils').discInvCdf2;
   rpois = require('./poisson').rpois;
   rgamma = require('./gamma').rgamma;

   /**
    * Evaluates the negative binomial distribution's density function at `x`,
    * where $x \geq 0$.
    *
    * `size` must be strictly positive, and `prob` must be in `(0, 1]`.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the
    * logarithm of the result.
    *
    * Adapted from dnbinom.c in R's code.
    *
    * @fullName dnbinom(size, prob, logp)(x)
    * @memberof nbinom
    */
   function dnbinom(size, prob, logp) {
      logp = logp === true;

      if (prob <= 0 || prob > 1 || size < 0) {
        return function(x) { return NaN; };
      }

      return function(x) {
         var p, ans;

         if (x < 0 || x !== Math.floor(x)) {
            return logp ? -Infinity : 0;
         }
         if (x === 0 && size === 0) {
            return logp ? 0 : 1;
         }

         ans = dbinomLog(x + size, prob)(size);
         p = size / (size + x);

         return logp ? Math.log(p) + ans : p * Math.exp(ans);
      };
   }

   /**
    * Evaluates the lower-tail cdf at `x` for the negative binomial distribution.
    *
    * `size` must be strictly positive, and `prob` must be in `(0, 1]`.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability instead.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    *
    * @fullName pnbinom(size, prob, lowerTail, logp)(x)
    * @memberof nbinom
    */
   function pnbinom(size, prob, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      return function(x) {

         if (isNaN(x + size + prob) || size < 0 || prob <= 0 || prob > 1) {
           return NaN;
         }
         if (size === 0) {
             return adjustLower(x >= 0 ? 1 : 0, lowerTail, logp);
         }
         if (x < 0) { return adjustLower(0, lowerTail, logp); }

         x = Math.floor(x + 1e-14);

         return beta.pbeta(size, x + 1, lowerTail, logp)(prob);
      };
   }

   /**
    * Evaluates the negative binomial distribution's quantile function
    * (inverse cdf) at `p`.
    *
    * `size` must be strictly positive, and `prob` must be in `(0, 1]`.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
    * interpreted as an upper tail probability (returns
    * $x$ such that $\textrm{prob}(X > x) = p)$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
    * the logarithm of the desired probability.
    *
    * @fullName qnbinom(size, prob, lowerTail, logp)(p)
    * @memberof nbinom
    */
   function qnbinom(size, prob, lowerTail, logp) {
      var mu, sigma, gamma;

      logp = logp === true;
      lowerTail = lowerTail !== false;

      mu = size * (1 - prob) / prob;
      sigma = Math.sqrt(size * (1 - prob) / prob / prob);
      gamma = (2 - prob) / prob / sigma;

      /* eslint-disable complexity */
      return pWrap(lowerTail, logp, function(ps) {
         var z, ret;

         if (isNaN(ps.p + prob + size)) { return NaN; }
         if (prob === 0 && size === 0) { return 0; }
         if (prob <= 0 || prob > 1 || size < 0) { return NaN; }
         if (prob === 1 || size === 0) { return 0; }

         if (ps.p === 0) { return 0; }
         if (ps.q === 0) { return Infinity; }

         z = ps.p > 0.9 ? qnorm(0, 1, false)(ps.q) : qnorm(0, 1)(ps.p); // initial value
         ret = Math.floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);

         // Supposed to do: ps.p *= 1 - 64*DBL_EPSILON here
         if (ps.p > 0.9) {
            return discInvCdf2(0, Infinity, ret, ps.q, pnbinom(size, prob, false));
         }
         return discInvCdf(0, Infinity, ret, ps.p, pnbinom(size, prob));
      });
      /* eslint-enable complexity */
   }

   /**
    * Returns a random variate from the negative binomial distribution.
    *
    * `size` must be strictly positive, and `prob` must be in `(0, 1]`.
    *
    * @fullName rnbinom(size, prob)()
    * @memberof nbinom
    */
   function rnbinom(size, prob) {
      var rg;

      if (size <= 0 || prob <= 0 || prob > 1) {
         return function() { return NaN; };
      }
      if (prob === 1) { return function() { return 0; }; }
      rg = rgamma(size, (1 - prob) / prob);
      return function() { return rpois(rg())(); };
   }

   return {
      /**
       * Returns an object representing a negative binomial distribution with
       * given parameters `size` and `prob`.  The object has
       * properties `d`, `p`, `q`, `r`.
       *
       * `size` must be strictly positive, and `prob` must be in `(0, 1]`.
       *
       * ```
       * nbinom(size, prob).d(x, logp)            // same as dnbinom(size, prob, logp)(x)
       * nbinom(size, prob).p(x, lowerTail, logp) // same as pnbinom(size, prob, lowerTail, logp)(x)
       * nbinom(size, prob).q(x, lowerTail, logp) // same as qnbinom(size, prob, lowerTail, logp)(x)
       * nbinom(size, prob).r()                   // same as rnbinom(size, prob)()
       * ```
       * @memberof nbinom
       */
      nbinom: function(size, prob) {
         return {
            d: function(x, logp) { return dnbinom(size, prob, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pnbinom(size, prob, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qnbinom(size, prob, lowerTail, logp)(p);
            },
            r: function() { return rnbinom(size, prob)(); }
         };
      },
      dnbinom: dnbinom,
      pnbinom: pnbinom,
      qnbinom: qnbinom,
      rnbinom: rnbinom
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
