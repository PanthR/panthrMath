(function(define) {'use strict';
define(function(require) {

   var C;

   C = {
      twopi: 2 * Math.PI,
      sqrt2pi: Math.sqrt(2 * Math.PI)
   };

   return Object.freeze(C);

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
