(function(define) {'use strict';
define(function(require) {

   /**
    * Implementation / polyfill for exp(x) - 1
    */

   return {
      expm1: Math.expm1 || function(x) { return Math.exp(x) - 1; }
   };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
