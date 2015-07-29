(function(define) {'use strict';
define(function(require) {

   /**
    * exp(x) - 1
    * Taken from:  Computation of the Incomplete Gamma Function Ratios
    * and their Inverse, by DiDonato and Morris.
    */
   var Rational, expm1;

   Rational = require('../rational');

   expm1 = (function() {
      var R = new Rational([
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
