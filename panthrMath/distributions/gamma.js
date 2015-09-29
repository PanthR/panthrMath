(function(define) {'use strict';
define(function(require) {

   // Gamma distribution

   var dpois;

   dpois = require('./poisson').dpois;

   // helper function
   function dgammaLog(a, s) {
      if (a < 1) {
         return function(x) {
            if (x === 0) { return Infinity; }
            return Math.log(a / x) + dpois(x / s, true)(a);
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
         return -Math.log(s) + dpois(x / s, true)(a - 1);
      };
   }

   // density function
   // 1 / (s^a * gamma(a)) * x ^ (a - 1) * exp(-x / s)
   function dgamma(a, s, logp) {
      logp = logp === true;
      return function(x) {
         return logp ? dgammaLog(a, s)(x) : Math.exp(dgammaLog(a, s)(x));
      };
   }

   return {
      dgamma: dgamma
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
