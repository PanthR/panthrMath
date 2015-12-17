(function(define) {'use strict';
define(function(require) {

   var C;

   C = {
      twopi: 2 * Math.PI,
      sqrt2pi: Math.sqrt(2 * Math.PI),
      log2pi: Math.log(2 * Math.PI),
      eulerGamma: 0.577215664901532860606512090082
   };

   return Object.freeze(C);

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
