(function(define) {'use strict';
define(function(require) {

   /**
    * Basic math functions.
    */

   var basicFunc;

   basicFunc = {};

   mixin(basicFunc,
      require('basicFunc/log1p'),
      require('basicFunc/expm1')
   );

   return basicFunc;

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
