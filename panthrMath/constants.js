(function(define) {'use strict';
define(function(require) {
   var C;

   /**
    * Constants module, provides frozen values for some frequently-used constants.
    * @module C
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   C = {
      /**
       * $2\pi$
       */
      twopi: 2 * Math.PI,
      /**
       *$\sqrt{2\pi}$
       */
      sqrt2pi: Math.sqrt(2 * Math.PI),
      /**
       * $\ln(2\pi)$
       */
      log2pi: Math.log(2 * Math.PI),
      /**
       * The Euler-Mascheroni constant
       * \\[\gamma = \lim \_ {n \to \infty} \left(-\ln(n) + \sum \_ {k=1}^n \frac{1}{k}\right) = 0.57721566 \ldots\\]
       *
       */
      eulerGamma: 0.577215664901532860606512090082
   };

   return Object.freeze(C);

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
