(function(define) {'use strict';
define(function(require) {

   var rgen;

   rgen = require('./rgen');

   // inverseCDF.js

   // takes a getMode function, updateLeft, and updateRight
   // The resulting function will call getMode with whatever arguments are
   // passed to it.  getMode returns the object with val = mode, prob =
   // corresponding probability.
   // updateLeft and updateRight will be called with `this` set to an object
   // with properties `val`, `prob`, and `update`.  They need to update the
   // object (for continuing the search, using the pdf for the distribution)
   // and return true iff the search can continue in that direction.
   return function(getMode, updateLeft, updateRight) {
      var mode, left, right;

      mode = getMode.apply(this, arguments);
      left = { update: updateLeft };
      right = { update: updateRight };

      return function() {
         var current, u;

         left.val = right.val = mode.val;
         left.prob = right.prob = mode.prob;
         current = right;

         u = rgen.random();
         do {
            if ( u <= current.prob ) { return current.val; }
            u -= current.prob;
            current = current === left ?
               getNext(right, left) :
               getNext(left, right);
         } while (u >= 0);
         throw new Error('total probability exceeds 1');
      };
   };

   function getNext(first, second) {
      if (first.update()) { return first; }
      if (second.update()) { return second; }
      throw new Error('both sides exhausted -- error!');
   }

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
