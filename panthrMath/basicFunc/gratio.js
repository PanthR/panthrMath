(function(define) {'use strict';
define(function(require) {

   /**
    * Implementation of the incomplete gamma function and its
    * complement.
    * See: "Computation of the Incomplete Gamma Function Ratios and their
    * Inverse", DiDonato and Morris
    *
    * P(a, x), Q(a, x) are the ratios, where P + Q = 1.
    * Domain: a >= 0, x >= 0, a + x != 0.
    */

   var C, phi, H, gratio, gratioc,
       Rational, Polynomial, lgamma, gamma, erf, erfc,
       Psmall, Qsmall, Pmedium, Qmedium, Plarge, Qlarge;

   Rational = require('../rational');
   Polynomial = require('../polynomial');
   lgamma = require('./lgamma').lgamma;
   gamma = require('./lgamma').gamma;
   erf = require('./erf').erf;
   erfc = require('./erf').erfc;
   C = require('../constants');

   /*
    * Computes x - 1 - ln(x)
    */
   phi = (function() {
      var phi1 = new Rational([
         0.00620886815375787,
        -0.224696413112536,
         0.333333333333333
      ], [
         0.354508718369557,
        -1.27408923933623,
         1
      ]);

      return function(x) {
         var r;
         if (x < 0) { return NaN; }
         if (x === 0) { return Infinity; }
         // TODO: Rethink this range based on paper's equation 6
         if (x < 0.82 || x > 1.18) { return x - 1 - Math.log(x); }

         r = (x - 1) / (x + 1);
         return 2 * r * r * ( 1 / (1 - r) - r * phi1.evalAt(r * r));
      };
   }());

   // Appendix C of DiDonato/Morris.
   // Computes 1/Gamma(x+1) - 1
   H = (function() {
      var w, w1;

      w = new Rational([
        -0.132674909766242e-3,
         0.266505979058923e-3,
         0.00223047661158249,
        -0.0118290993445146,
         0.930357293360349e-3,
         0.118378989872749,
        -0.244757765222226,
        -0.771330383816272,
        -0.422784335098468
      ], [
         0.0559398236957378,
         0.273076135303597,
         1
      ]);
      w1 = new Rational([
         0.589597428611429e-3,
        -0.00514889771323592,
         0.00766968181649490,
         0.0597275330452234,
        -0.230975380857675,
        -0.409078193005776,
         0.577215664901533
      ], [
         0.00423244297896961,
         0.0261132021441447,
         0.158451672430138,
         0.427569613095214,
         1
      ]);

      return function(x) {
         if (x <= -1) { return NaN; }
         if (x < -0.5 || x > 1.5) { return 1 / gamma(x + 1) - 1; }
         if (x <= 0) { return x * (1 + w.evalAt(x)); }
         if (x <= 0.5) { return x * w1.evalAt(x); }
         if (x <= 1) { return (x - 1) / x * w.evalAt(x - 1); }
         return (x - 1) / x * (w1.evalAt(x - 1) - 1);
      };
   }());

   // Psmall = (function() {

   //    return function(a) {
   //       return function(x) {

   //       };
   //    };
   // }());

   // Curried form of gamma ratio  (P(a, x))
   function gratio(a) {
      if (a == 0.5) { return function(x) { return erf(Math.sqrt(x)); } }
      if (a < 1) { return Psmall(a); }
      if (a < 20) { return Pmedium(a); }
      Plarge(a);
   }
   function gratioc(a) {
      if (a == 0.5) { return function(x) { return erfc(Math.sqrt(x)); } }
      if (a < 1) { return Qsmall(a); }
      if (a < 20) { return Qmedium(a); }
      Qlarge(a);
   }

   return {
      gratio: gratio,
      gratioc: gratioc
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
