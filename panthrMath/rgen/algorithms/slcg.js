(function(define) {'use strict';
define(function(require) {

   function slcg(a, c, m) {
      var x, q, r;             // x is the seed
      q = Math.floor(m / a);
      r = m % a;
      function random() {
         // take care to ensure exact computations, depending on value of a
         // x = (a * x + c) % m might cause overflow; instead...
         var k;
         k = Math.floor(x / q);
         x = a * (x - k * q) - k * r + c;
         while (x < 0) { x = x + m; }  // this is an 'if' in the book?
         return x / m;
      }
      function setSeed(i) {
         x = i;
         return this;
      }
      return {
         random: random,
         setSeed: setSeed
      };
   }

   return slcg;
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
