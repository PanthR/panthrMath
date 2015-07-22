(function(define) {'use strict';
define(function(require) {

   /**
    * Implementations / polyfills for basic math functions.
    */

   return {
      log1p: Math.log1p || function(x) { return Math.log(1 + x); },
      expm1: Math.expm1 || function(x) { return Math.exp(x) - 1; }
   };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
