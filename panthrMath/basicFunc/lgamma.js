(function(define) {'use strict';
define(function(require) {

   var C, lgamma, lgammaLanczos, lgammaNear1or2, Rational;

   Rational = require('../rational');

   C = require('../constants');

   /* The Lanczos version works well for positive inputs not close to 1. */
   lgammaLanczos = (function() {
      // z Must be positive
      var cs = [
         0.99999999999980993227684700473478,
         676.520368121885098567009190444019,
         -1259.13921672240287047156078755283,
         771.3234287776530788486528258894,
         -176.61502916214059906584551354,
         12.507343278686904814458936853,
         -0.13857109526572011689554707,
         9.984369578019570859563e-6,
         1.50563273514931155834e-7
      ];
      return function(z) {
         var t, i, ser;

         t = z + 7.5;
         t -= (z + 0.5) * Math.log(t);
         ser = cs[0];
         for (i = 1; i < cs.length; i += 1) {
            ser += cs[i] / (z + i);
         }

         return -t + Math.log(C.sqrt2pi * ser / z);
      };
   }());

   /**
    * lgamma(x) for 0.8 <= x <= 2.25
    * Taken from:  Computation of the Incomplete Gamma Function Ratios
    * and their Inverse, by DiDonato and Morris.
    */
   lgammaNear1or2 = (function() {
      var R1, R2;
      R1 = new Rational([
         -0.00271935708322958,
         -0.0673562214325671,
         -0.402055799310489,
         -0.780427615533591,
         -0.168860593646662,
          0.844203922187225,
          0.577215664901533
         ], [
          0.667465618796164e-3,
          0.0325038868253937,
          0.361951990101499,
          1.56875193295039,
          3.12755088914843,
          2.88743195473681,
          1
         ]);
      R2 = new Rational([
          0.497958207639485e-3,
          0.0170502484022650,
          0.156513060486551,
          0.565221050691933,
          0.848044614534529,
          0.422784335098467
         ], [
          0.116165475989616e-3,
          0.00713309612391000,
          0.101552187439830,
          0.548042109832463,
          1.24313399877507,
          1
         ]);
      return function(x) {
         if (x <= 1.6) { return (1 - x) * R1.evalAt(x - 1); }
         return (x - 2) * R2.evalAt(x - 2);
      };

   }());

   lgamma = function(x) {
      if (x === 1 || x === 2) { return 0; }
      if (x >= 0.8 && x <= 2.25) { return lgammaNear1or2(x); }
      // todo handle negative case
      return lgammaLanczos(x);
   };

   return { lgamma: lgamma };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
