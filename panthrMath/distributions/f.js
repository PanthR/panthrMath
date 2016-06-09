(function(define) {
'use strict';
define(function(require) {

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the Fisher-Snedecor F distribution with parameters `df1` and `df2`
    * (degrees of freedom).  It is the distribution of the ratio of the mean
    * squares of `df1` and `df2` independent standard normal variates.
    *
    * `df` provides access to the probability density function,
    * `pf` to the cumulative distribution function, `qf` to the
    * quantile function (inverse cdf)
    * and `rf` to random deviates.
    *
    * Finally, you can use `fdistr` to obtain an object
    * representing the F distribution for given values of
    * `df1` and `df2`.
    *
    * @module distributions.f
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var gamma, beta, dbinomLog, chisq, utils;

   dbinomLog = require('../basicFunc/dbinomLog').dbinomLog;
   gamma = require('./gamma');
   beta = require('./beta');
   chisq = require('./chisq');
   utils = require('../utils');

   /**
    * Evaluates the F distribution's density function at `x`, where
    * $x \geq 0$.
    *
    * Expects `df1` and `df2` to be positive.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the
    * logarithm of the result.
    *
    * Adapted from df.c in R's code.
    *
    * @fullName df(df1, df2, logp)(x)
    * @memberof f
    */
   function df(df1, df2, logp) {
      logp = logp === true;

      if (df1 <= 0 || df2 <= 0 || utils.hasNaN(df1, df2)) {
        return function(x) { return NaN; };
      }

      /* eslint-disable complexity */
      return function(x) {
         var f, q, p, dens;

         if (utils.hasNaN(x)) { return NaN; }
         if (x < 0) { return logp ? -Infinity : 0; }
         if (x === 0) {
            if (df1 > 2) { return logp ? -Infinity : 0; }
            if (df1 === 2) { return logp ? 0 : 1; }
            return Infinity;
         }
         if (df1 === Infinity && df2 === Infinity) {
            return x === 1 ? Infinity : logp ? -Infinity : 0;
         }
         if (df2 === Infinity) { return gamma.dgamma(df1 / 2, 2 / df1, logp)(x); }
         if (df1 > 1e14) {
            return logp ? gamma.dgamma(df2 / 2, 2 / df2, logp)(1 / x) - 2 * Math.log(x)
                        : gamma.dgamma(df2 / 2, 2 / df2, logp)(1 / x) / (x * x);
         }

         f = 1 / (df2 + x * df1);
         q = df2 * f;
         p = x * df1 * f;

         if (df1 >= 2) {
            p = x === Infinity ? 1 : p;
            dens = dbinomLog((df1 + df2 - 2) / 2, p)((df1 - 2) / 2);
            f = df1 * q / 2;
         } else {
            dens = dbinomLog((df1 + df2) / 2, p)(df1 / 2);
            f = df1 * df1 * q / (2 * p * (df1 + df2));
         }
         return logp ? Math.log(f) + dens : f * Math.exp(dens);
      };
      /* eslint-enable complexity */
   }

   /**
    * Evaluates the lower-tail cdf at `x` for the F distribution.
    *
    * Expects `df1` and `df2` to be positive.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability instead.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    *
    * @fullName pf(df1, df2, lowerTail, logp)(x)
    * @memberof f
    */
   function pf(df1, df2, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (df1 <= 0 || df2 <= 0 || utils.hasNaN(df1, df2)) {
        return function(x) { return NaN; };
      }

      return function(x) {
         if (utils.hasNaN(x)) { return NaN; }
         if (x <= 0) { return utils.adjustLower(0, lowerTail, logp); }

         // deal with df === INF -- use pchisqr
         if (df2 === Infinity) {
            if (df1 === Infinity) {
               return x < 1 ? utils.adjustLower(0, lowerTail, logp)
                  : x === 1 ? logp ? -Math.log(2) : 0.5
                            : utils.adjustLower(1, lowerTail, logp);
            }
            return chisq.pchisq(df1, lowerTail, logp)(x * df1);
         }
         if (df1 === Infinity) {
            return chisq.pchisq(df2, !lowerTail, logp)(df2 / x);
         }

         if (df1 * x > df2) {
            return beta.pbeta(df2 / 2, df1 / 2,
                              !lowerTail, logp)(df2 / (df2 + df1 * x));
         }

         return beta.pbeta(df1 / 2, df2 / 2,
                           lowerTail, logp)(df1 * x / (df2 + df1 * x));
      };
   }

   /**
    * Evaluates the F distribution's quantile function
    * (inverse cdf) at `p`.
    *
    * Expects `df1` and `df2` to be positive.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
    * interpreted as an upper tail probability (returns
    * $x$ such that $\textrm{prob}(X > x) = p)$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
    * the logarithm of the desired probability.
    *
    * @fullName qf(df1, df2, lowerTail, logp)(p)
    * @memberof f
    */
   function qf(df1, df2, lowerTail, logp) {
      var qbeta;

      logp = logp === true;
      lowerTail = lowerTail !== false;
      qbeta = beta.qbeta(df2 / 2, df1 / 2, !lowerTail, logp);

      if (df1 <= 0 || df2 <= 0 || utils.hasNaN(df1, df2)) {
         return function(x) { return NaN; };
      }

      return utils.qhelper(lowerTail, logp, 0, Infinity, function(p) {
         if (df1 <= df2 && df2 > 4e5) {
            return df1 === Infinity ? 1
                   : chisq.qchisq(df1, lowerTail, logp)(p) / df1;
         }
         if (df1 > 4e5) { return df2 / chisq.qchisq(df2, !lowerTail, logp)(p); }
         return (1 / qbeta(p) - 1) * (df2 / df1);
      });
   }

   /**
    * Returns a random variate from the F distribution, computed as a ratio
    * of chi square variates.
    *
    * Expects `df1` and `df2` to be positive.
    *
    * @fullName rf(df1, df2)()
    * @memberof f
    */
   function rf(df1, df2) {
      var rchisqRatio;

      rchisqRatio = function(d) { return chisq.rchisq(d)() / d; };

      if (df1 <= 0 || df2 <= 0) { return function() { return NaN; }; }

      return function() {
         return rchisqRatio(df1) / rchisqRatio(df2);
      };
   }

   return {
      /**
       * Returns an object representing an F distribution with
       * given degrees of freedom `df1` and `df2`.  The object has
       * properties `d`, `p`, `q`, `r`.
       * ```
       * fdistr(df1, df2).d(x, logp)            // same as df(df1, df2, logp)(x)
       * fdistr(df1, df2).p(x, lowerTail, logp) // same as pf(df1, df2, lowerTail, logp)(x)
       * fdistr(df1, df2).q(x, lowerTail, logp) // same as qf(df1, df2, lowerTail, logp)(x)
       * fdistr(df1, df2).r()                   // same as rf(df1, df2)()
       * ```
       * @memberof f
       */
      fdistr: function(df1, df2) {
         return {
            d: function(x, logp) { return df(df1, df2, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pf(df1, df2, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qf(df1, df2, lowerTail, logp)(p);
            },
            r: function() { return rf(df1, df2)(); }
         };
      },
      df: df,
      pf: pf,
      qf: qf,
      rf: rf
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
