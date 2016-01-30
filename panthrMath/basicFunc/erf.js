(function(define) {
'use strict';
define(function(require) {

   /*

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

   /**
    * Implementation for erf, the real error function:
    * $$\textrm{erf}(x) = \frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2} dt$$
    * Based on: *Rational Chebyshev Approximations for the Error
    * Function*, by W. J. Cody, 1969
    * @memberof basicFunc
    */
   erf = function(x) {
      if (x < -0.5) { return -erf(-x); }
      if (x < 0.5) { return x * R4small.evalAt(x * x); }
      return 1 - erfc(x);
   };

   /**
    * Implementation for erfc, the complement of the real error function:
    * $$\textrm{erfc}(x) = \frac{2}{\sqrt{\pi}}\int_x^\infty e^{-t^2} dt$$
    * Based on: *Rational Chebyshev Approximations for the Error
    * Function*, by W. J. Cody, 1969
    * @memberof basicFunc
    */
   erfc = function(x) {
      var x2inv;

      if (x < -0.5) { return 2 - erfc(-x); }
      if (x < 0.5) { return 1 - erf(x); }
      if (x <= 4) { return Math.exp(-x * x) * R8med.evalAt(x); }
      x2inv = 1 / (x * x);
      return Math.exp(-x * x) / x *
         (piSqrtInv + x2inv * R5large.evalAt(x2inv));
   };

   // erf__ = (function(){
   //    var Rat1, Rat2, Rat3, c;
   //    c = 0.564189583547756;

   //    Rat1 = Rational.new([
   //       7.7105849500132e-5,
   //      -.00133733772997339,
   //       .0323076579225834,
   //       .0479137145607681,
   //       .128379167095513 + 1
   //    ], [
   //       .00301048631703895,
   //       .0538971687740286,
   //       .375795757275549,
   //       1
   //    ]);

   //    Rat2 = Rational.new([
   //       -1.36864857382717e-7,
   //       .564195517478974,
   //       7.21175825088309,
   //       43.1622272220567,
   //       152.98928504694,
   //       339.320816734344,
   //       451.918953711873,
   //       300.459261020162
   //    ], [
   //       1.,
   //       12.7827273196294,
   //       77.0001529352295,
   //       277.585444743988,
   //       638.980264465631,
   //       931.35409485061,
   //       790.950925327898,
   //       300.459260956983
   //    ]);

   //    Rat3 = Rational.new([
   //       2.10144126479064,
   //       26.2370141675169,
   //       21.3688200555087,
   //       4.6580782871847,
   //       .282094791773523
   //    ], [
   //       94.153775055546,
   //       187.11481179959,
   //       99.0191814623914,
   //       18.0124575948747,
   //       1
   //    ]);

   //    return function(x) {
   //       var ax, sgnx;

   //       ax = Math.abs(x);
   //       sgnx = x < 0 ? -1 : 1;
   //       if (ax <= 0.5) { return x * Rat1.evalAt(x * x); }
   //       if (ax <= 4) {
   //          return (1 - Math.exp(-x * x) * Rat2.evalAt(ax)) * sgnx;
   //       }
   //       if (ax >= 5.8) { return sgnx; }
   //       return (1 - Math.exp(-x * x) *
   //                   (c - Rat3.evalAt(1 / x / x) / x / x) / ax
   //              ) * sgnx;
   //    };
   // }());

   return {
      erf: erf,
      erfc: erfc
   };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
