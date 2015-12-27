(function(define) {'use strict';
define(function(require) {

   /**
    * Probability Distributions
    *
    * TODO
    *
    * @module distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var distributions, mixin;

   mixin = require('./utils.js').mixin;

   distributions = mixin({},
      /** Constructor for finite distributions. See `module:finite`. */
      require('./distributions/finite'),
      /** The normal distribution. See `module:normal`. */
      require('./distributions/normal'),
      /** The gamma distribution. See `module:gamma`. */
      require('./distributions/gamma'),
      /** The Poisson distribution. See `module:poisson`. */
      require('./distributions/poisson'),
      /** The t distribution. See `module:t`. */
      require('./distributions/t'),
      /** The uniform distribution. See `module:uniform`. */
      require('./distributions/uniform'),
      /** The beta distribution. See `module:beta`. */
      require('./distributions/beta'),
      /** The binomial distribution. See `module:binomial`. */
      require('./distributions/binom')
   );

   return distributions;

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
