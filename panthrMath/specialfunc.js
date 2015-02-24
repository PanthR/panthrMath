(function(define) {'use strict';
define(function(require) {

   var C, lgamma, lfactorial, lchoose;

   C = require('./constants');

   lgamma = (function() {
      // z Must be positive
      var cs = [
         0.99999999999980993227684700473478,
         676.520368121885098567009190444019, -1259.13921672240287047156078755283,
         771.3234287776530788486528258894, -176.61502916214059906584551354,
         12.507343278686904814458936853, -0.13857109526572011689554707,
         9.984369578019570859563e-6, 1.50563273514931155834e-7
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

   // Expects x >= 0;
   lfactorial = function lfactorial(x) {
      return x === 0 ? 0 : lgamma(x + 1);
   };

   lchoose = function lchoose(n, k) {
      return lfactorial(n) - lfactorial(k) - lfactorial(n - k);
   };

   return {
      lgamma: lgamma,
      lfactorial: lfactorial,
      lchoose: lchoose
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
