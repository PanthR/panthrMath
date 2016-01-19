(function(define) {'use strict';
define(function(require) {

   /**
    * Random number generation
    *
    * panthrMath implements its own random number generation, to provide
    * consistency across platforms.  Currently, only one algorithm is provided.
    * See: *Random Number Generation and Monte Carlo Methods*, 2nd edition,
    * by J. Gentle (Springer).
    *
    * Available algorithms:
    * - `slcg`: simple linear congruential generator, using the
    *     congruence $x\_i \equiv a x _ {i-1} + c \bmod m$ where
    *     $a = 16807$, $c = 0$, and $m = 2^{31} - 1$.
    *
    * @module rgen
    * @version 0.0.1
    * @author Haris Skiadas <skiadas@hanover.edu>
    * Barb Wahl <wahl@hanover.edu>
    */

   /*
    Some thoughts...
    Each algorithm is for a uniform distribution on [0, 1].  Different
    algorithms exist with various tradeoffs.

    an algorithm needs to be able to provide:
       - settings for use at creation
       - setSeed function
       - random function

    rgen will be a single object with
       - algorithms property
       - setAlgorithm function
       - distributions property (uniform, normal, etc.) -- for example,
       rgen.distributions.normal takes mean and standard deviation and returns
       a random number from the corresponding normal distribution
       - setSeed function (from the current algorithm)
       - random function (from the current algorithm)

    once a specific algorithm has been set, commands such as
    rgen.setSeed or rgen.random will use that algorithm
   */

   var rgen, slcg;

   rgen = {
      algorithms: {},
      /**
       * Sets a seed for the random number generator.  Expects `seed` to
       * be an integer.
       */
      setSeed: function(seed) {
         throw new Error('need to select an algorithm first');
         // return rgen;
      },
      /**
       * Generates a random number uniformly on the interval $(0, 1)$.
       */
       random: function() {
         throw new Error('need to select an algorithm first');
      },
      /**
       * Sets the algorithm to be used for random number generation.
       * Expects `name` to be one of the strings returned by `getAlgorithms`.
       */
      setAlgorithm: function(name) {
         // "name" encapsulates some simple options like precision
         rgen.random = rgen.algorithms[name].random;
         rgen.setSeed = rgen.algorithms[name].setSeed;
         return rgen.setRandomSeed();
      },
      /**
       * Returns a list of the available algorithms for random number
       * generation.
       */
      getAlgorithms: function() {
         return Object.keys(rgen.algorithms);
      }
   };

   /**
    * Sets a seed for random number generation, based on the system date.
    * @fullName setRandomSeed()
    */
   rgen.setRandomSeed = function() {
      return rgen.setSeed(new Date());
   };

   // Load Algorithms
   slcg = require('./algorithms/slcg');
   rgen.algorithms.slcg = slcg(16807, 0, Math.pow(2, 31) - 1);
   rgen.setAlgorithm('slcg');

   return rgen;
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
