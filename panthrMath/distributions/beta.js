(function(define) {'use strict';
define(function(require) {

   var bratio;

   bratio = require('../basicFunc/bratio').bratio;

   function pbeta(a, b) {
      return function(x) {
         if (a <= 0 || b <= 0) { return NaN; }
         if (x <= 0) { return 0; }
         if (x >= 1) { return 1; }
         return bratio(a, b, x);
      };
   }

   pbeta.log = function(a, b) {
      return function(x) {
         if (a <= 0 || b <= 0) { return NaN; }
         if (x <= 0) { return -Infinity; }
         if (x >= 1) { return 0; }
         return bratio.log(a, b, x);
      };
   };

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
