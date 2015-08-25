(function(define) {'use strict';
define(function(require) {

   // Gamma distribution

   var C, twopi, sqrt2pi, dpoisLog;

   C = require('../constants');
   twopi = C.twopi;
   sqrt2pi = C.sqrt2pi;
   dpoisLog = require('./poisson').dpois.log;

   // helper function
   function dgammaLog(a, s) {
      if (a < 1) {
         return function(x) {
            if (x === 0) { return Infinity; }
            return Math.log(a / x) + dpoisLog(x / s)(a);
         };
      }
      if (a === 1) {
         return function(x) {
            return -x / s - Math.log(s);
         };
      }
      // a > 1
      return function(x) {
         if (x === 0) { return -Infinity; }
         return -Math.log(s) + dpoisLog(x / s)(a - 1);
      };
   }

   // density function
   // 1 / (s^a * gamma(a)) * x ^ (a - 1) * exp(-x / s)
   function dgamma(a, s) {
      return function(x) {
         return dgammaLog(a, s)(x);
      };
   }

   dgamma.log = dgammaLog;

   return {
      dgamma: dgamma
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
