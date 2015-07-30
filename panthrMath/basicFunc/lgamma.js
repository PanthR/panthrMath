(function(define) {'use strict';
define(function(require) {

   /**
    * Implementation of the logarithm of the gamma function.
    * Only works for positive x values.
    */

   var C, lgamma, lgammaLanczos, lgammaNear1or2, Rational, Polynomial;

   Rational = require('../rational');
   Polynomial = require('../polynomial');

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

   // Performs (2,2)-pade approximation. (See lgamma.c in gsl)
   // Corr is a Polynomial
   function pade22(n1, n2, d1, d2, c, corr) {
      return function(eps) {
         var eps5, pade;
         pade = c * (eps + n1) * (eps + n2) / ((eps + d1) * (eps + d2));
         eps5 = eps * eps * eps * eps * eps;
         return eps * (pade + eps5 * corr.evalAt(eps));
      };
   }

   /**
    * lgamma(x) for 0.8 <= x <= 2.25
    * Taken from:  lgamma.c (GSL)
    */
   lgammaNear1or2 = (function() {
      var p1, p2;
      p1 = pade22(
         -1.0017419282349508699871138440, // n1
          1.7364839209922879823280541733, // n2
          1.2433006018858751556055436011, // d1
          5.0456274100274010152489597514, // d2
          2.0816265188662692474880210318, // c
         new Polynomial([
             0.03141928755021455,
            -0.02594027398725020,
             0.01931961413960498,
            -0.01192457083645441,
             0.004785324257581753
         ])
      );
      p2 = pade22(
         1.000895834786669227164446568, // n1
         4.209376735287755081642901277, // n2
         2.618851904903217274682578255, // d1
         10.85766559900983515322922936, // d2
         2.85337998765781918463568869,  // c
         new Polynomial([
          0.0000407220927867950,
         -0.0000693271800931282,
          0.0001067287169183665,
         -0.0001365435269792533,
          0.0001139406357036744
         ])
      );

      return function(x) {
         return x <= 1.6 ? p1(x - 1) : p2 (x - 2);
      };

   }());

   lgamma = function(x) {
      if (x <= 0) { return NaN; }
      if (x === 1 || x === 2) { return 0; }
      if (x >= 0.8 && x <= 2.25) { return lgammaNear1or2(x); }
      return lgammaLanczos(x);
   };

   return { lgamma: lgamma };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
