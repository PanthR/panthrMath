(function(define) {'use strict';
define(function(require) {

   /*
    * Math computation library for PanthR
    * @module panthrMath
    * @version 0.0.1
    * @author Haris Skiadas <skiadas@hanover.edu>
    * Barb Wahl <wahl@hanover.edu>
    */

   var panthrMath;

   function mixin(target) {
      Array.prototype.slice.call(arguments, 1)
         .forEach(function(source) {
            Object.keys(source).forEach(function(key) {
               target[key] = source[key];
            });
         });
      return target;
   }

   // Functions

   panthrMath = {};
   panthrMath.C = require('./panthrMath/constants');

   mixin(panthrMath,
      require('./panthrMath/specialfunc')
   );

   return panthrMath;
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
