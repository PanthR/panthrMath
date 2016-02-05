(function(define) {
'use strict';
define(function(require) {

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the chi square distribution ... TODO
    *
    * @module distributions.chisq
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var gamma;

   gamma = require('./gamma');

   /**
    * TODO
    *
    * @fullName dchisq(df, logp)(x)
    * @memberof chisq
    */
   function dchisq(df, logp) {
      return gamma.dgamma(df / 2, 2, logp);
   }

   /**
    * TODO
    *
    * @fullName pchisq(df, lowerTail, logp)(x)
    * @memberof chisq
    */
   function pchisq(df, lowerTail, logp) {
      return gamma.pgamma(df / 2, 2, lowerTail, logp);
   }

   /**
    * TODO
    *
    * @fullName qchisq(df, lowerTail, logp)(p)
    * @memberof chisq
    */
   function qchisq(df, lowerTail, logp) {
      return gamma.qgamma(df / 2, 2, lowerTail, logp);
   }

   /**
    * TODO
    *
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
