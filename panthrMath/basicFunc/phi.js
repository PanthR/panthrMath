(function(define) {'use strict';
define(function(require) {

   /*
    * Appendix C of DiDonato/Morris.
    * Computes x - 1 - ln(x)
    */

   var phi1, Rational;

   Rational = require('../rational');

   phi1 = new Rational([
      0.00620886815375787,
     -0.224696413112536,
      0.333333333333333
   ], [
      0.354508718369557,
     -1.27408923933623,
      1
   ]);

   function phi(x) {
      var rat;
      if (x < 0) { return NaN; }
      if (x === 0) { return Infinity; }
      if (x < 0.82 || x > 1.18) { return x - 1 - Math.log(x); }

      rat = (x - 1) / (x + 1);
      return 2 * rat * rat * ( 1 / (1 - rat) - rat * phi1.evalAt(rat * rat));
   }

   return { phi: phi };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
