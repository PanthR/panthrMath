(function(define) {'use strict';
define(function(require) {

   return {
      mixin: function mixin(target) {
         Array.prototype.slice.call(arguments, 1)
            .forEach(function(source) {
               Object.keys(source).forEach(function(key) {
                  target[key] = source[key];
               });
            });
         return target;
      }
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
