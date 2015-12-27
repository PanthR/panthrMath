(function(define) {'use strict';
define(function(require) {

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the continuous uniform distribution.
    *
    * @module distributions.uniform
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var rgen;

   rgen = require('../rgen/rgen');

   /**
    * TODO
    * @memberof uniform
    */
   function dunif(min, max, logp) {
      logp = logp === true;

      if (!(min < max)) { return function(x) { return NaN; }; }

      return function(x) {
         var p;
         p = x < min || x > max ? 0 : 1 / (max - min);

         return logp ? Math.log(p) : p;
      };
   }

   /**
    * TODO
    * @fullName punif(min, max, lowerTail, logp)(x)
    * @memberof uniform
    */
   function punif(min, max, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (!(min < max)) { return function(x) { return NaN; }; }

      return function(x) {
         var p;

         p = x <= min ? 0 :
             x >= max ? 1 :
                        (x - min) / (max - min);
         p = lowerTail ? p : 1 - p;

         return logp ? Math.log(p) : p;
      };
   }

   /**
    * TODO
    * @fullName qunif(min, max, lowerTail, logp)(p)
    * @memberof uniform
    */
   function qunif(min, max, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (!(min < max)) { return function(x) { return NaN; }; }

      return function(p) {
         p = logp ? Math.exp(p) : p;

         if (p < 0 || p > 1) { return NaN; }
         return lowerTail ? min + p * (max - min) : max - p * (max - min);
      };
   }

   /**
    * TODO
    * @memberof uniform
    */
   function runif(min, max) {
      return function() {
         return min + (max - min) * rgen.random();
      };
   }

   return {
      /**
       * Returns an object representing a uniform distribution, with properties `d`, `p`, `q`, `r`.
       * ```
       * unif(a, b).d(x, logp)            // same as dunif(a, b, logp)(x)
       * unif(a, b).p(x, lowerTail, logp) // same as punif(a, b, lowerTail, logp)(x)
       * unif(a, b).q(x, lowerTail, logp) // same as qunif(a, b, lowerTail, logp)(x)
       * unif(a, b).r(n)                  // same as runif(a, b)(n)
       * ```
       * @memberof uniform
       */
      unif: function(min, max) {
         return {
            d: function(x, logp) { return dunif(min, max, logp)(x); },
            p: function(q, lowerTail, logp) {
               return punif(min, max, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qunif(min, max, lowerTail, logp)(p);
            },
            r: function(n) { return runif(min, max)(n); }
         };
      },
      dunif: dunif,
      punif: punif,
      qunif: qunif,
      runif: runif
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
