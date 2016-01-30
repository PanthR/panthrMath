(function(define) {
'use strict';
define(function(require) {

   var Rational, log1p;

   Rational = require('../rational');

   /**
    * Computes $$\textrm{log1p}(x) = \ln(1+x)$$ with reasonable precision near $x = 0$.
    *
    * Based on:  *Computation of the Incomplete Gamma Function Ratios
    * and their Inverse*, by DiDonato and Morris, 1986
    * @fullName log1p(x)
    * @memberof basicFunc
    */
   log1p = (function() {
      var R;

      R = Rational.new([
         -0.0178874546012214,
          0.405303492862024,
         -1.29418923021993,
          1
         ], [
         -0.0845104217945565,
          0.747811014037616,
         -1.62752256355323,
          1
         ]);
      return function(x) {
         var t;

         if (Math.abs(x) > 0.375) { return Math.log(1 + x); }
         t = x / (x + 2);
         return 2 * t * R.evalAt(t * t);
      };
   }());

   return { log1p: log1p };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
