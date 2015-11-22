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
