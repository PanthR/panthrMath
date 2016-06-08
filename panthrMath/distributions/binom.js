(function(define) {
'use strict';
define(function(require) {

   /**
    * Provides probability mass function, cumulative distribution function,
    * quantile function, and random number generator
    * for the Binomial distribution, which is defined by the pmf
    * $$f(x;n,p) = {n \choose x}p^x(1-p)^{n-x}$$
    * where $x=0,\ldots,n$, $n$ is the *size* and $p$ is the probability of success.
    *
    * `dbinom` provides access to this probability mass function,
    * `pbinom` to the cumulative distribution function, `qbinom` to the
    * quantile function (inverse cdf)
    * and `rbinom` to random deviates.
    *
    * Finally, you can use `binom` to obtain an object
    * representing the distribution for some values of the parameters.
    * @module distributions.binomial
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var dbinomLog, pbeta, qnorm, pWrap, discInvCdf, discInvCdf2, inverseCDF, utils;

   dbinomLog = require('../basicFunc/dbinomLog').dbinomLog;
   pbeta = require('./beta').pbeta;
   qnorm = require('./normal').qnorm;
   pWrap = require('../utils').pWrap;
   utils = require('../utils');
   discInvCdf = require('../utils').discInvCdf;
   discInvCdf2 = require('../utils').discInvCdf2;
   inverseCDF = require('../rgen/inverseCDF');

   /**
    * Returns the Binomial probability at `x`:
    * $$\textrm{dbinom}(\textrm{size}, p)(x) = \binom{\textrm{size}}{p}p^{x}(1-p)^{(\textrm{size}-x)}$$
    * where `x` is an integer, $0 \leq x \leq \textrm{size}$.
    *
    * `size` is a positive integer (number of trials) and $0 \leq p \leq 1$
    * (the probability of success on a single trial).
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the
    * logarithm of the result.
    *
    * Based on: *Fast and Accurate Computation of Binomial Probabilities*,
    * by Catherine Loader, 2000
    *
    * @fullName dbinom(size, p, logp)(x)
    * @memberof binomial
    */
   function dbinom(size, p, logp) {
      var lbinom;

      logp = logp === true;
      lbinom = dbinomLog(size, p);

      if (utils.hasNaN(size, p) ||
          size < 0 || p < 0 || p > 1 ||
          size !== Math.floor(size)) {
         return function() { return NaN; };
      }

      return function(x) {
         if (utils.hasNaN(x)) { return NaN; }

         if (x < 0 || x > size || Math.floor(x) !== x ||
                      utils.isInfinite(x)) {
            return logp ? -Infinity : 0;
         }

         return logp ? lbinom(x) : Math.exp(lbinom(x));
      };
   }

   /**
    * Evaluates the Binomial cumulative distribution
    * function at `x` (lower tail probability):
    * $$\textrm{pbinom}(\textrm{size}, p)(x) = \sum_{k \leq x} \binom{\textrm{size}}{p}p^{k}(1-p)^{(\textrm{size}-k)}$$
    *
    * `size` is a positive integer (number of trials) and $0 \leq p \leq 1$
    * (the probability of success on a single trial).
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability instead:
    * $$\textrm{pbinom}(\textrm{size}, p, \textrm{false})(x) = \sum_{k>x} \binom{\textrm{size}}{p}p^{k}(1-p)^{(\textrm{size}-k)}$$
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    *
    * @fullName pbinom(size, p, lowerTail, logp)(x)
    * @memberof binomial
    */
   function pbinom(size, p, lowerTail, logp) {
      var goodParams;

      logp = logp === true;
      lowerTail = lowerTail !== false;

      goodParams = size >= 0 && size < Infinity && size === Math.floor(size) &&
                   p >= 0 && p <= 1;
      if (!goodParams) { return function(x) { return NaN; }; }
      return function(x) {
         var res;

         if (isNaN(x)) { return NaN; }
         res = lowerTail === x < 0 ? -Infinity : 0;
         x = Math.floor(x);
         if (x >= 0 && x < size) {
            res = pbeta(x + 1, size - x, !lowerTail, true)(p);
         }

         return logp ? res : Math.exp(res);
      };
   }
   /**
    * Returns the quantile corresponding to `prob`
    * for the $\textrm{Binomial}(size, p)$ distribution.
    * In general, for a discrete probability
    * distribution, the *quantile* is defined as the smallest domain value
    * `x` such that $F(x) \geq prob$, where $F$ is the cumulative
    * distribution function.
    *
    * `size` is a positive integer (number of trials) and $0 \leq p \leq 1$
    * (the probability of success on a single trial).
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `prob` is
    * interpreted as an upper tail probability.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `prob` as
    * the logarithm of the desired probability.
    *
    * @fullName qbinom(size, prob, lowerTail, logp)(prob)
    * @memberof binomial
    */
   function qbinom(size, prob, lowerTail, logp) {
      var goodParams, mu, sigma, gamma;

      logp = logp === true;
      lowerTail = lowerTail !== false;
      goodParams = size >= 0 && size < Infinity && size === Math.floor(size) &&
                   prob >= 0 && prob <= 1;

      if (!goodParams) { return function(p) { return NaN; }; }
      if (prob === 0) {
         return pWrap(lowerTail, logp,
            function(ps) { return ps.q === 0 ? size : 0; }
         );
      }
      if (prob === 1) {
         return pWrap(lowerTail, logp,
            function(ps) { return ps.p === 0 ? 0 : size; }
         );
      }
      // 0 < prob < 1
      mu = size * prob;
      sigma = Math.sqrt(mu * (1 - prob));
      gamma = (1 - 2 * prob) / sigma;

      return utils.qhelper(lowerTail, logp, 0, size, function(p) {
         var z, ret;

         if (!utils.isFinite(size) || !utils.isFinite(prob)) { return NaN; }
         if (size === 0 || prob === 0) { return 0; }

         z = qnorm(0, 1, lowerTail, logp)(p); // initial value
         ret = Math.floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);
         if (ret > size) { ret = size; }

         if (lowerTail) {
            return discInvCdf(0, size, ret, p, pbinom(size, prob, true, logp));
         }

         return discInvCdf2(0, size, ret, p, pbinom(size, prob, false, logp));
      });
   }

   // rbinom, using inverseCDF

   /**
    * Returns a random variate from the $\textrm{Binomial}(size, p)$ distribution.
    *
    * `size` is a positive integer (number of trials) and $0 \leq p \leq 1$
    * (the probability of success on a single trial).
    *
    * @fullName rbinom(size, p)()
    * @memberof binomial
    */
   function rbinom(size, p) {
      var mode;

      mode = { val: Math.floor(size * p) };
      mode.prob = dbinom(size, p)(mode.val);

      return inverseCDF(
         function getMode() { return mode; },
         function updateLeft() {
            if (this.val === 0) { return false; }
            this.prob *= this.val / (size - this.val + 1) * ((1 - p) / p);
            this.val -= 1;
            return true;
         },
         function updateRight() {
            if (this.val === size) { return false; }
            this.prob *= (size - this.val) / (this.val + 1) * (p / (1 - p));
            this.val += 1;
            return true;
         }
      );
   }

   return {
      /**
       * Returns an object representing a binomial distribution, with properties `d`, `p`, `q`, `r`.
       * ```
       * binom(size, p).d(x, logp)            // same as dbinom(size, p, logp)(x)
       * binom(size, p).p(x, lowerTail, logp) // same as pbinom(size, p, lowerTail, logp)(x)
       * binom(size, p).q(x, lowerTail, logp) // same as qbinom(size, p, lowerTail, logp)(x)
       * binom(size, p).r()                   // same as rbinom(size, p)()
       * ```
       * @memberof binomial
       */
      binom: function(size, p) {
         return {
            d: function(x, logp) { return dbinom(size, p, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pbinom(size, p, lowerTail, logp)(q);
            },
            q: function(prob, lowerTail, logp) {
               return qbinom(size, p, lowerTail, logp)(prob);
            },
            r: function() { return rbinom(size, p)(); }
         };
      },
      dbinom: dbinom,
      pbinom: pbinom,
      qbinom: qbinom,
      rbinom: rbinom
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
