(function(define) {
'use strict';
define(function(require) {

   var Rational, expm1;

   Rational = require('../rational');

   /**
    * Computes $$\textrm{expm1}(x) = e^x - 1$$ with reasonable precision near $x = 0$.
    *
    * Based on:  *Computation of the Incomplete Gamma Function Ratios
    * and their Inverse*, by DiDonato and Morris, 1986
    * @fullName expm1(x)
    * @memberof basicFunc
    */
   expm1 = (function() {
      var R;

      R = Rational.new([
          0.0238082361044469,
          0.914041914819518e-9,
          1
         ], [
          0.595130811860248e-3,
         -0.0119041179760821,
          0.107141568980644,
         -0.499999999085958,
          1
         ]);
      return function(x) {
         if (x < -0.15) { return Math.exp(x) - 1; }
         if (x > 0.15) { return Math.exp(x) * (1 - 1 / Math.exp(x)); }
         return x * R.evalAt(x);
      };
   }());

   return { expm1: expm1 };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
