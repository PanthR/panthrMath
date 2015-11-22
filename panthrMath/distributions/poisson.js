(function(define) {'use strict';
define(function(require) {

   // Poisson distribution
   // No input validation provided.

   var lpoisson, pgamma;

   lpoisson = require('../basicFunc/lpoisson').lpoisson;
   pgamma = require('./gamma').pgamma;

   // density / pdf
   function dpois(lambda, logp) {
      logp = logp === true;
      return function(x) {
         return logp ? lpoisson(lambda)(x) : Math.exp(lpoisson(lambda)(x));
      };
   }

   function ppois(lambda, lowerTail, logp) {
      lowerTail = lowerTail !== false;
      logp = logp === true;

      if (!(lambda >= 0)) { return function(x) { return NaN; }; }

      return function(x) {
         var ret;

         if (x >= 0 && lambda > 0) {
            return pgamma(Math.floor(x + 1e-10) + 1, 1, !lowerTail, logp)(lambda);
         }

         ret = lowerTail ? 1 : 0;
         if (x < 0) { ret = lowerTail ? 0 : 1; }
         return logp ? Math.log(ret) : ret;
      };
   }

   return {
      pois: function(lambda) {
         return {
            d: function(x, logp) { return dpois(lambda, logp)(x); },
            p: function(q, lowerTail, logp) {
               return ppois(lambda, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qpois(lambda, lowerTail, logp)(p);
            },
            r: function(n) { return rpois(lambda)(n); }
         };
      },
      dpois: dpois,
      ppois: ppois
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
