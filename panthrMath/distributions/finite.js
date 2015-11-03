(function(define) {'use strict';
define(function(require) {

   var utils;

   utils = require('../utils');

   /**
    * Finite distribution
    * The object o must have either
    * - properties `xs`, `ws`, which are arrays of equal length, or
    * - properties `f`, `min`, `max`, where the available values are the
    * sequential numbers from min to max with a step of 1, and f(i) gives
    * the probability of the value i.
    *
    * The xs are assumed to be distinct and in increasing order.
    *
    * The ws are treated as "weights": They need to be nonnegative, and
    * if they do not add up to 1 they will be rescaled.
    */
   function finite(o) {
      var sum, xs, ws, cums;
      if (typeof o.f === 'function') {
         o = populateObjectFromFunction(o);
      }
      xs = o.xs;
      ws = o.ws;

      sum = ws.reduce(function(v, acc) {
         if (!(v >= 0)) { throw new Error('finite: Negative weights'); }
         return v + acc;
      }, 0);
      if (!utils.relativelyCloseTo(sum, 1)) {
         ws = ws.map(function(v) { return v / sum; });
      }
      sum = 0;  // Reusing sum for the cumulative probabilities (cdf)
      cums = ws.map(function(v) { sum += v; return sum; });

      return {
         d: function(x, logp) {
            var i, ret;

            logp = logp === true;
            if (x < xs[0] || x > xs[xs.length - 1]) {
               ret = 0;
            } else {
               i = binSearch(xs, function(y) { return y <= x; });
               ret = xs[i] === x ? ws[i] : 0;
            }
            return logp ? Math.log(ret) : ret;
         },
         p: function(q, lowerTail, logp) {
            var i, ret;
            logp = logp === true;
            lowerTail = lowerTail !== false;
            if (isNaN(q)) { return NaN; }
            if (q < xs[0]) {
               ret = 0;
            } else if (q >= xs[xs.length - 1]) {
               ret = 1;
            } else {
               i = binSearch(xs, function(y) { return y <= q; });
               ret = cums[i];
            }
            ret = lowerTail ? ret : 1 - ret;
            return logp ? Math.log(ret) : ret;
         },
         q: function(p, lowerTail, logp) {
            logp = logp === true;
            lowerTail = lowerTail !== false;
            return utils.pWrap(lowerTail, logp, function(prob) {
               if (isNaN(prob)) { return NaN; }
               if (prob < 0 || prob > 1) { return NaN; }
               return xs[binSearch2(cums, function(y) { return y >= prob; })];
            })(p);
         },
         r: function(n) {
            // TODO: When we bring rgen in
         }
      };
   }

   /*
    * Given a strictly increasing array of xs, and a function pred that is a predicate
    * and assuming that there is an x in the set of xs such that:
    * pred(y) is true iff y <= x
    * this function then finds and returns the index in the array of x
    */
   function binSearch(xs, pred) {
      var a, b, mid;
      a = 0;
      b = xs.length;
      while (b - a > 1) {
         mid = Math.floor((a + b) / 2);
         if (pred(xs[mid])) {
            a = mid;
         } else {
            b = mid;
         }
      }
      return a;
   }
   /*
    * Given a strictly increasing array of xs, and a function pred that is a predicate
    * and assuming that there is an x in the set of xs such that:
    * pred(y) is true iff y >= x
    * this function then finds and returns the index in the array of x
    */
   function binSearch2(xs, pred) {
      var a, b, mid;
      a = -1;
      b = xs.length - 1;
      while (b - a > 1) {
         mid = Math.floor((a + b) / 2);
         if (pred(xs[mid])) {
            b = mid;
         } else {
            a = mid;
         }
      }
      return b;

   }

   function populateObjectFromFunction(o) {
      var i, obj;

      obj = { xs: [], ws: [] };
      if (!(o.min <= o.max)) {
         throw new Error('finite: Invalid range specification, expected min <= max');
      }
      for (i = o.min; i <= o.max; i += 1) {
         obj.xs.push(i);
         obj.ws.push(o.f(i));
      }

      return obj;
   }

   return {
      finite: finite,
      dfinite: function dfinite(o, logp) {
         var distr = finite(o);
         return function(x) { return distr.d(x, logp); };
      },
      pfinite: function pfinite(o, lowerTail, logp) {
         var distr = finite(o);
         return function(q) { return distr.p(q, lowerTail, logp); };
      },
      qfinite: function qfinite(o, lowerTail, logp) {
         var distr = finite(o);
         return function(p) { return distr.q(p, lowerTail, logp); };
      },
      rfinite: function rfinite(o) {
         return finite(o).r;
      }
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
