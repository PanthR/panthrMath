(function(define) {'use strict';
define(function(require) {

   var Polynomial = require('./polynomial');

   // creates a rational function
   // num can be a Polynomial or an array of coefficients for the poly constr.
   // same for denom
   function Rational(num, denom) {
      this.num = new Polynomial(num);
      this.denom = new Polynomial(denom);
   }

   Rational.prototype.evalAt = function evalAt(x) {
      return this.num.evalAt(x) / this.denom.evalAt(x);
   };


   return Rational;

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
