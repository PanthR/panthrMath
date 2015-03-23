(function(define) {'use strict';
define(function(require) {

   var twopi = require('../constants').twopi;

   // log density
   function dnormLog(mu, sigma) {
      var c = Math.log(sigma * sigma * twopi);
      return function(x) {
         var z = (x - mu) / sigma;
         return -0.5 * (c + z * z);
      };
   }

   // density / pdf
   function dnorm(mu, sigma) {
      var f = dnormLog(mu, sigma);
      return function(x) { return Math.exp(f(x)); };
   }

   // cdf
   function pnorm(mu, sigma) {
      return function(x) {

      };
   }

   // inverse cdf
   function qnorm(mu, sigma) {
      return function(p) {

      };
   }

   return {
      dnormLog: dnormLog,
      dnorm: dnorm,
      pnorm: pnorm,
      qnorm: qnorm
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
