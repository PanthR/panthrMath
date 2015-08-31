(function(define) {'use strict';
define(function(require) {

   /**
    * bd0
    * Taken from:  Loader 2000
    */
   function bd0(x, np) {
      var ej, j, s, sTemp, v;
      if (Math.abs(x - np) < 0.1 * (x + np)) {
         v = (x - np) / (x + np);
         s = (x - np) * v;
         ej = 2 * x * v;
         for (j = 1; j < 10000; j += 1) {
            ej *= v * v;
            sTemp = s + ej / (2 * j + 1);
            if (sTemp === s) {
               return sTemp;
            }
            s = sTemp;
         }
         throw new Error('bd0(x, np) made it to 10000 ...');
      }
      return x * Math.log(x / np) + np - x;
   }

   return { bd0: bd0 };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
