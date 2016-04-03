(function(define) {
'use strict';
define(function(require) {

   /**
    * TODO: Update the documentation
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the normal distribution, which is defined by the pdf
    * $$f(x;\mu,\sigma) = \frac{1}{\sigma \sqrt{2\pi}} e^{\displaystyle -\frac{(x-\mu)^2}{2\sigma^2}}$$
    * where $x\in(-\infty,\infty)$ and as usual $\mu$ is the mean and $\sigma>0$
    * is the standard deviation.
    *
    * `dnorm` provides access to this probability density function,
    * `pnorm` to the cumulative distribution function, `qnorm` to the
    * quantile function (inverse cdf)
    * and `rnorm` to random deviates.
    *
    * Finally, you can use `normal` to obtain an object
    * representing the distribution for some values of the parameters.
    * @module distributions.normal
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var normal, C, logRoot2pi, recipRoot2pi;

   normal = require('./normal');
   C = require('./constants');
   logRoot2pi = C.log2pi * 0.5;
   recipRoot2pi = 1 / C.sqrt2pi;

   /**
    * TODO
    * assuming sdlog > 0
    * support for dlnorm is x  > 0
    * @fullName dnorm(meanlog, sdlog, logp)(x)
    * @memberof normal
    */
   function dlnorm(meanlog, sdlog, logp) {
      logp = logp === true;
      return function(x) {
         var z;

         if (x <= 0) { return logp ? -Infinity : 0; }

         z = (Math.log(x) - meanlog) / sdlog;
         return logp ? -(logRoot2pi + 0.5 * z * z + Math.log(x * sdlog))
                     : recipRoot2pi * Math.exp(-0.5 * z * z) / (x * sdlog);
      };
   }

   /**
    * TODO FIX ME
    * Evaluates the lower-tail cdf at `x` for the Normal distribution:
    * $$\textrm{pnorm}(\mu, \sigma)(x) = \frac{1}{\sigma \sqrt{2\pi}}\int_{-\infty}^x e^{\displaystyle -\frac{(t-\mu)^2}{2\sigma^2}}dt$$
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
    * the upper tail probability instead.
    *
    * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
    * of the result.
    *
    * Expects $\sigma > 0$.
    *
    * @fullName pnorm(meanlog, sdlog, lowerTail, logp)(x)
    * @memberof normal
    */
   function plnorm(meanlog, sdlog, lowerTail, logp) {
      var pnorm;

      lowerTail = lowerTail !== false;
      logp = logp === true;
      pnorm = normal.pnorm(meanlog, sdlog, lowerTail, logp);

      return function(x) {
         if (x <= 0) {
            return lowerTail ? logp ? -Infinity : 0
                             : logp ? 0 : 1;
         }

         return pnorm(Math.log(x));
      };
   }

   /**
    * TODO
    * Evaluates the Normal distribution's quantile function (inverse cdf) at `p`:
    * $$\textrm{qnorm}(\mu, \sigma)(p) = x \textrm{ such that } \textrm{prob}(X \leq x) = p$$
    * where $X$ is a random variable with the $N(\mu,\sigma)$ distribution.
    *
    * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
    * interpreted as an upper tail probability (returns
    * $x$ such that $\textrm{prob}(X > x) = p)$.
    *
    * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
    * the logarithm of the desired probability.
    *
    * Expects $\sigma > 0$.
    * @fullName qnorm(meanlog, sdlog, lowerTail, logp)(p)
    * @memberof normal
    */
   function qlnorm(meanlog, sdlog, lowerTail, logp) {
      var qnorm;

      logp = logp === true;
      lowerTail = lowerTail !== false;
      qnorm = normal.qnorm(meanlog, sdlog, lowerTail, logp);

      if (sdlog <= 0) { return function(x) { return NaN; }; }

      return function(p) {
         return Math.exp(qnorm(p));
      };
   }

   /**
    * TODO
    * Returns a random variate from the $N(\mu, \sigma)$ distribution.
    *
    * Expects $\sigma > 0$.
    *
    * Uses a rejection polar method.
    * @fullName rnorm(meanlog, sdlog)()
    * @memberof normal
    */
   function rlnorm(meanlog, sdlog) {
      var rnorm;

      rnorm = normal.rnorm(meanlog, sdlog);

      return function() {
         return Math.exp(rnorm());
      };
   }

   return {
      /**
       * Returns an object representing a normal distribution, with properties `d`, `p`, `q`, `r`.
       * ```
       * lognormal(meanlog, sdlog).d(x, logp)            // same as dlnorm(meanlog, sdlog, logp)(x)
       * lognormal(meanlog, sdlog).p(x, lowerTail, logp) // same as plnorm(meanlog, sdlog, lowerTail, logp)(x)
       * lognormal(meanlog, sdlog).q(x, lowerTail, logp) // same as qlnorm(meanlog, sdlog, lowerTail, logp)(x)
       * lognormal(meanlog, sdlog).r()                   // same as rlnorm(meanlog, sdlog)()
       * ```
       * @memberof lognormal
       */
      lognormal: function(meanlog, sdlog) {
         return {
            d: function(x, logp) { return dlnorm(meanlog, sdlog, logp)(x); },
            p: function(q, lowerTail, logp) {
               return plnorm(meanlog, sdlog, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qlnorm(meanlog, sdlog, lowerTail, logp)(p);
            },
            r: function() { return rlnorm(meanlog, sdlog)(); }
         };
      },
      dlnorm: dlnorm,
      plnorm: plnorm,
      qlnorm: qlnorm,
      rlnorm: rlnorm
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
