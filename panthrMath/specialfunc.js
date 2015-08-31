(function(define) {'use strict';
define(function(require) {

   var lgamma, lfactorial, lchoose;

   lgamma = require('./basicFunc/lgamma').lgamma;

   // Expects x >= 0;
   lfactorial = function lfactorial(x) {
      return x === 0 ? 0 : lgamma(x + 1);
   };

   lchoose = function lchoose(n, k) {
      return lfactorial(n) - lfactorial(k) - lfactorial(n - k);
   };

   return {
      lfactorial: lfactorial,
      lchoose: lchoose
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
