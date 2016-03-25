(function(define) {
'use strict';
define(function(require) {

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the chi squared distribution on $k$ degrees of freedom, defined by
    * the pdf:
    * $$f(x; u) = \frac{1}{2^u\Gamma(u)}x^{u - 1} e^{-x/2}$$
    * where $u = k / 2$, the integer $k > 0$ is the degrees of freedom
    * and $x \in [0, \infty)$.
    *
    * `dchisq` provides access to this probability density function,
    * `pchisq` to the cumulative distribution function, `qchisq` to the
    * quantile function (inverse cdf)
    * and `rchisq` to random deviates.
    *
    * Finally, you can use `chisq` to obtain an object
    * representing the chi squared distribution for a given value of
    * the degrees of freedom $k$.
    *
    * The implementation uses the gamma distribution. See `module:gamma`.
    * @module distributions.chisq
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var gamma;

   gamma = require('./gamma');

   /**
    * Evaluates the chi squared distribution's density function at `x`.
    *
    * `df` is the degrees of freedom $k$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the
    * logarithm of the result.
    *
    * @fullName dchisq(df, logp)(x)
    * @memberof chisq
    */
   function dchisq(df, logp) {
      return gamma.dgamma(df / 2, 2, logp);
   }

   /**
    * Evaluates the lower-tail cdf at `x` for the chi squared distribution:
    * $$\textrm{pchisq}(k)(x) = \frac{\gamma(k/2, x/2)}{\Gamma(k/2)}$$
    *
    * where $\gamma$ is the incomplete gamma function, and
    * `df` is the degrees of freedom $k$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability instead.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    *
    * @fullName pchisq(df, lowerTail, logp)(x)
    * @memberof chisq
    */
   function pchisq(df, lowerTail, logp) {
      return gamma.pgamma(df / 2, 2, lowerTail, logp);
   }

   /**
    * Evaluates the chi-squared distribution's quantile function
    * (inverse cdf) at `p`:
    * $$\textrm{qchisq}(k)(p) = x \textrm{ such that } \textrm{prob}(X \leq x) = p$$
    * where $X$ is a random variable with the chi-squared distribution.
    *
    * `df` is the degrees of freedom $k$.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
    * interpreted as an upper tail probability (returns
    * $x$ such that $\textrm{prob}(X > x) = p)$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
    * the logarithm of the desired probability.
    *
    * @fullName qchisq(df, lowerTail, logp)(p)
    * @memberof chisq
    */
   function qchisq(df, lowerTail, logp) {
      return gamma.qgamma(df / 2, 2, lowerTail, logp);
   }

   /**
    * Returns a random variate from the chi squared distribution with `df` degrees
    * of freedom.
    *
    * @fullName rchisq(df)()
    * @memberof chisq
    */
   function rchisq(df) {
      return gamma.rgamma(df / 2, 2);
   }

   return {
      /**
       * Returns an object representing a chi squared distribution for `df` degrees
       * of freedom, with properties `d`, `p`, `q`, `r`.
       * ```
       * chisq(df).d(x, logp)            // same as dchisq(df, logp)(x)
       * chisq(df).p(x, lowerTail, logp) // same as pchisq(df, lowerTail, logp)(x)
       * chisq(df).q(x, lowerTail, logp) // same as qchisq(df, lowerTail, logp)(x)
       * chisq(df).r()                   // same as rchisq(df)()
       * ```
       * @memberof chisq
       */
      chisq: function(df) {
         return {
            d: function(x, logp) { return dchisq(df, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pchisq(df, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qchisq(df, lowerTail, logp)(p);
            },
            r: function() { return rchisq(df)(); }
         };
      },
      dchisq: dchisq,
      pchisq: pchisq,
      qchisq: qchisq,
      rchisq: rchisq
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
