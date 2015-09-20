(function(define) {'use strict';
define(function(require) {

   /**
    * Implementation for erf and erfc, the real error function and its
    * complement, using Rational Chebyshev Approximations for the Error
    * Function by W. J. Cody.
    */

   var erf, erfc, R4small, R8med, R5large, Rational, piSqrtInv;

   piSqrtInv = 1 / Math.sqrt(Math.PI);
   Rational = require('../rational');

   /* R4small is Cody's function R_lm with l = m = 4 for |x| <= 0.5*/
   R4small = Rational.new([
         0.1857777061846031526730,
         3.161123743870565596947,
         113.8641541510501556495,
         377.4852376853020208137,
         3209.377589138469472562
      ], [
         1,
         23.60129095234412093499,
         244.0246379344441733056,
         1282.616526077372275645,
         2844.236833439170622273
      ]);

   /* R8med is Cody's function R_lm with l = m = 8 for
      0.5 <= x <= 4.0
   */
   R8med = Rational.new([
         2.15311535474403846343e-8,
         0.564188496988670089180,
         8.88314979438837594118,
         66.1191906371416294775,
         298.635138197400131132,
         881.952221241769090411,
         1712.04761263407058314,
         2051.07837782607146532,
         1230.33935479799725272
      ], [
         1,
         15.7449261107098347253,
         117.693950891312499305,
         537.181101862009857509,
         1621.38957456669018874,
         3290.79923573345962678,
         4362.61909014324715820,
         3439.36767414372163696,
         1230.33935480374942043
      ]);

   /* R5large is Cody's function R_lm with l = m = 5 for x >= 4 */
   R5large = Rational.new([
         -0.0163153871373020978498,
         -0.305326634961232344035,
         -0.360344899949804439429,
         -0.125781726111229246204,
         -0.0160837851487422766278,
         -6.58749161529837803157e-4
      ], [
         1,
         2.56852019228982242072,
         1.87295284992346047209,
         0.527905102951428412248,
         0.0605183413124413191178,
         0.00233520497626869185443
      ]);

   erf = function erf(x) {
      if (x < -0.5) { return -erf(-x); }
      if (x < 0.5) { return x * R4small.evalAt(x * x); }
      return 1 - erfc(x);
   };

   erfc = function erfc(x) {
      var x2inv;
      if (x < -0.5) { return 2 - erfc(-x); }
      if (x < 0.5) { return 1 - erf(x); }
      if (x <= 4) { return Math.exp(-x * x) * R8med.evalAt(x); }
      x2inv = 1 / (x * x);
      return Math.exp(-x * x) / x *
         (piSqrtInv + x2inv * R5large.evalAt(x2inv));
   };

   return {
      erf: erf,
      erfc: erfc
   };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
