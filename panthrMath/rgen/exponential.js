(function(define) {'use strict';
define(function(require) {

   var rgen;

   rgen = require('./rgen');

   // exponential.js

   // Uses inverse cdf
   return function(lambda) {
      return function() {
         return -Math.log(rgen.random()) / lambda;
      };
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
