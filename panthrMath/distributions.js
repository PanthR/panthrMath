(function(define) {'use strict';
define(function(require) {

   var distributions, mixin;

   mixin = require('./utils.js').mixin;

   /**
    * Probability Distributions
    * @module distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   distributions = mixin({},
      /**
       * TODO
       */
      require('./distributions/finite'),
      /**
       * TODO
       */
      require('./distributions/normal'),
      /**
       * TODO
       */
      require('./distributions/gamma'),
      /**
       * TODO
       */
      require('./distributions/poisson'),
      /**
       * TODO
       */
      require('./distributions/t'),
      /**
       * TODO
       */
      require('./distributions/uniform'),
      /**
       * TODO
       */
      require('./distributions/beta'),
      /**
       * TODO
       */
      require('./distributions/binom')
   );

   return distributions;

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
