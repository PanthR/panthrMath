(function(define) {'use strict';
define(function(require) {

   // Appendix C of DiDonato/Morris.
   // Computes 1/Gamma(x+1) - 1

   var w, w1, gamma, Rational;

   gamma = require('./lgamma').gamma;
   Rational = require('../rational');

   w = new Rational([
     -0.132674909766242e-3,
      0.266505979058923e-3,
      0.00223047661158249,
     -0.0118290993445146,
      0.930357293360349e-3,
      0.118378989872749,
     -0.244757765222226,
     -0.771330383816272,
     -0.422784335098468
   ], [
      0.0559398236957378,
      0.273076135303597,
      1
   ]);
   w1 = new Rational([
      0.589597428611429e-3,
     -0.00514889771323592,
      0.00766968181649490,
      0.0597275330452234,
     -0.230975380857675,
     -0.409078193005776,
      0.577215664901533
   ], [
      0.00423244297896961,
      0.0261132021441447,
      0.158451672430138,
      0.427569613095214,
      1
   ]);

   function gam1(x) {
      if (x <= -1) { return NaN; }
      if (x < -0.5 || x > 1.5) { return 1 / gamma(x + 1) - 1; }
      if (x <= 0) { return x * (1 + w.evalAt(x)); }
      if (x <= 0.5) { return x * w1.evalAt(x); }
      if (x <= 1) { return (x - 1) / x * w.evalAt(x - 1); }
      return (x - 1) / x * (w1.evalAt(x - 1) - 1);
   }

   return { gam1: gam1 };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
