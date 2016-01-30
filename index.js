(function(define) {
'use strict';
define(function(require) {

   /**
    * Math computation library for PanthR
    *
    * Provides special functions, probability distributions,
    * and random number generators; inspired by the R code.
    * @module panthrMath
    * @version 0.0.1
    * @noPrefix
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */

   var panthrMath, mixin, rgen;

   mixin = require('./panthrMath/utils').mixin;

   // Functions
   panthrMath = {};
   /** Frequently used constants.  See `module:C`
    */
   panthrMath.C = require('./panthrMath/constants');
   /* Abstraction for polynomial functions. */
   panthrMath.Polynomial = require('./panthrMath/polynomial');
   /* Abstraction for rational functions. */
   panthrMath.Rational = require('./panthrMath/rational');

   rgen = require('./panthrMath/rgen/rgen');
   ['setAlgorithm', 'getAlgorithms',
    'setSeed', 'setRandomSeed', 'random'
   ].forEach(function(key) { panthrMath[key] = rgen[key]; });

   mixin(panthrMath,
      /** A collection of basic and special functions. See `module:basicFunc`. */
      require('./panthrMath/basicfunc'),
      /** Various probability distributions. See `module:distributions`. */
      require('./panthrMath/distributions')
   );

   return panthrMath;
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
