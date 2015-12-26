(function(define) {'use strict';
define(function(require) {

   /**
    * Provides density function, constructor, cdf, inverse cdf, and random
    * number generator for finite discrete probability distributions.
    *
    * @module distributions.finite
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var utils, rgen;

   utils = require('../utils');
   rgen = require('../rgen/rgen');

   /*
    * Finite distribution
    * The object o must have either
    * - properties `xs`, `ws`, which are arrays of equal length, or
    * - properties `f`, `min`, `max`, where the available values are the
    * sequential numbers from min to max with a step of 1, and f(i) gives
    * the probability of the value i.
    *
    * The xs are assumed to be distinct and in increasing order.
    *
    * The ws are treated as "weights": They need to be _positive_, and
    * if they do not add up to 1 they will be rescaled.
    *
    * q tries to invert p but cannot be an exact inverse.  In particular,
    * if asked for the quantile <=  which the left area is 0, q will return
    * min and if asked for the quantile <= which the left area is 1,
    * q returns max.  The edge cases for !lowerTail are symmetrical to these:
    * q will return 0 or max for right-tail area of 1 or 0, respectively.
    */
    /**
     * TODO
     * @memberof finite
     */
   function finite(o) {
      var sum, xs, ws, cumsLeft, cumsRight, i, finObj;
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
      cumsLeft = ws.map(function(v) { sum += v; return sum; });
      cumsRight = [];
      cumsRight[ws.length - 1] = 0;
      for (i = ws.length - 2; i >= 0; i -= 1) {
         cumsRight[i] = cumsRight[i + 1] + ws[i + 1];
      }
      finObj = {
         d: function(x, logp) {
            var ind, ret;

            logp = logp === true;
            if (x < xs[0] || x > xs[xs.length - 1]) {
               ret = 0;
            } else {
               ind = binSearch(xs, function(y) { return y <= x; });
               ret = xs[ind] === x ? ws[ind] : 0;
            }
            return logp ? Math.log(ret) : ret;
         },
         /* eslint-disable complexity */
         p: function(q, lowerTail, logp) {
            var ind, ret;
            logp = logp === true;
            lowerTail = lowerTail !== false;
            if (isNaN(q)) { return NaN; }
            if (q < xs[0]) {
               ret = lowerTail ? 0 : 1;
            } else if (q >= xs[xs.length - 1]) {
               ret = lowerTail ? 1 : 0;
            } else {
               ind = binSearch(xs, function(y) { return y <= q; });
               ret = lowerTail ? cumsLeft[ind] : cumsRight[ind];
            }
            ret = Math.min(ret, 1);
            return logp ? Math.log(ret) : ret;
         },
         q: function(p, lowerTail, logp) {
            logp = logp === true;
            lowerTail = lowerTail !== false;
            if (logp) { p = Math.exp(p); }
            if (isNaN(p)) { return NaN; }
            if (p < 0 || p > 1) { return NaN; }
            if (p === 1) { return lowerTail ? xs[xs.length - 1] : xs[0]; }
            if (p === 0) { return lowerTail ? xs[0] : xs[xs.length - 1]; }
            if (lowerTail) {
               // we need the smallest x such that the area left and equal
               // to x is >= p
               return xs[binSearch2(cumsLeft, function(y) {
                  return utils.relativelyCloseTo(y, p, 1e-14) || y > p;
               })];
            }
            // !lowerTail
            // we need the smallest x such that the area strictly above
            // x is <= p
            return xs[binSearch2(cumsRight, function(y) {
               return utils.relativelyCloseTo(y, p, 1e-14) || y <= p;
            })];
         },
         /* eslint-enable complexity */
         r: function() {
            return finObj.q(rgen.random());
         }
      };

      return finObj;
   }


   /*
    * Given an array of xs, and a predicate on those values,
    * supposing the predicate is true on some initial portion of the
    * array and false elsewhere, the function returns the largest index
    * where the predicate is true.  If the predicate is always false,
    * returns 0.
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
    * Given an array of xs, and a predicate on those values,
    * supposing the predicate is false on some initial portion of the
    * array and true elsewhere, the function returns the smallest index
    * where the predicate is true.  If the predicate is always false,
    * returns the last index.
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
      dfinite:
      /**
       * TODO
       * @memberof finite
       */
       function dfinite(o, logp) {
         var distr = finite(o);
         return function(x) { return distr.d(x, logp); };
      },
      pfinite:
      /**
       * TODO
       * @memberof finite
       */
       function pfinite(o, lowerTail, logp) {
         var distr = finite(o);
         return function(q) { return distr.p(q, lowerTail, logp); };
      },
      qfinite:
      /**
       * TODO
       * @memberof finite
       */
       function qfinite(o, lowerTail, logp) {
         var distr = finite(o);
         return function(p) { return distr.q(p, lowerTail, logp); };
      },
      rfinite:
      /**
       * TODO
       * @memberof finite
       */
       function rfinite(o) {
         return finite(o).r;
      }
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
