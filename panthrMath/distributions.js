(function(define) {'use strict';
define(function(require) {

   // distributions.js

   var distributions, mixin;

   mixin = require('./utils.js').mixin;

   distributions = mixin({},
      require('./distributions/finite'),
      require('./distributions/normal'),
      require('./distributions/gamma'),
      require('./distributions/poisson'),
      require('./distributions/t'),
      require('./distributions/uniform'),
      require('./distributions/beta'),
      require('./distributions/binom')
   );

   return distributions;

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
