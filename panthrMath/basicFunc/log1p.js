(function(define) {'use strict';
define(function(require) {

   /**
    * Implementation / polyfill for log(1 + x)
    */

   return {
      log1p: Math.log1p || function(x) { return Math.log(1 + x); }
   };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
