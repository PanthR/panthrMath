(function(define) {'use strict';
define(function(require) {

   /* repeat, for carrying out repeated improvement until stopping
    * condition is satisfied.
    *
    * `init` is the initial value,
    * `step` is a function which returns the new value, and
    * `n` is a non-negative integer (number of steps to take),
    *  a function returning a boolean, or missing.
    */
   function repeat(init, step, n) {
      var prev, curr, done;
      curr = init;
      if (n === 0) { return curr; }
      done = typeof n === 'function' ? n :
             n > 0 ? function() { n -= 1; return n === 0; } :
                     function() { return curr === prev; };
      while (!done()) {
         prev = curr;
         curr = step();
      }
      return curr;
   }

   /* mixin */
   return {
      mixin: function mixin(target) {
         Array.prototype.slice.call(arguments, 1)
            .forEach(function(source) {
               Object.keys(source).forEach(function(key) {
                  target[key] = source[key];
               });
            });
         return target;
      },
      repeat: repeat
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
