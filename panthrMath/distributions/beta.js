(function(define) {'use strict';
define(function(require) {

   var bratio, dbeta, qbeta;

   bratio = require('../basicFunc/bratio').bratio;

   function pbeta(a, b, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      return function(x) {
         var lp;
         if (a <= 0 || b <= 0) { return NaN; }
         if (x > 0 && x < 1) {
            return bratio(a, b, x, lowerTail, logp);
         }
         lp = x <= 0 ? 0 : 1;
         lp = lowerTail ? lp : 1 - lp;
         return logp ? Math.log(lp) : lp;
      };
   }

   return {
      dbeta: dbeta,
      pbeta: pbeta,
      qbeta: qbeta
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
