(function(define) {'use strict';
define(function(require) {
   var C;

   /**
    * Constants module
    * @module C
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   C = {
      /**
       * TODO
       */
      twopi: 2 * Math.PI,
      /**
       * TODO
       */
      sqrt2pi: Math.sqrt(2 * Math.PI),
      /**
       * TODO
       */
      log2pi: Math.log(2 * Math.PI),
      /**
       * TODO
       */
      eulerGamma: 0.577215664901532860606512090082
   };

   return Object.freeze(C);

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
