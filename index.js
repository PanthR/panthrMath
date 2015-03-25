(function(define) {'use strict';
define(function(require) {

   /*
    * Math computation library for PanthR
    * @module panthrMath
    * @version 0.0.1
    * @author Haris Skiadas <skiadas@hanover.edu>
    * Barb Wahl <wahl@hanover.edu>
    */

   var panthrMath, mixin;

   mixin = require('./panthrMath/utils').mixin;

   // Functions

   panthrMath = {};
   panthrMath.C = require('./panthrMath/constants');
   panthrMath.Polynomial = require('./panthrMath/polynomial');
   panthrMath.Rational = require('./panthrMath/rational');

   mixin(panthrMath,
      require('./panthrMath/specialfunc'),
      require('./panthrMath/distributions')
   );

   return panthrMath;
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
