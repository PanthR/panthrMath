(function(define) {'use strict';
define(function(require) {

   var basicfunc;

   basicfunc = require('../basicfunc');

   /*
    * bratio -- see Algorithm 708
    *
    * Options object `opt` can have properties `lowerTail`, `log`:
    * `lowerTail` defaults to true, `log` defaults to false.
    */
    function bratio(a, b, x, opt) {

      // By convention, if 0 < x < 1 then a == 0 -> result = 1
      // and b == 0 -> result = 0. Can't have a == b == 0.
    }



   return {
      bratio: bratio
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
