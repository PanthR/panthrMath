(function(define) {'use strict';
define(function(require) {

   // distributions.js

   var distributions, mixin;

   mixin = require('./utils.js').mixin;

   distributions = mixin({},
      {}
   );

   return distributions;

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
