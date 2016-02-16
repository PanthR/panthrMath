(function(define) {
'use strict';
define(function(require) {

   /**
    * Provides probability mass function, constructor, cumulative distribution
    * function, quantile function, and random number generator for finite
    * discrete probability distributions.
    *
    * For all members of the `module:finite` module, distributions are provided
    * via appropriate properties of the object `o`. There are two possible options
    * for that object:
    * - It can have properties `xs`, `ws`, which are arrays of equal length,
    * corresponding to the values the distribution can take and their
    * respective probabilities. In this case the `xs` are assumed to be distinct and
    * in increasing order, and the `ws` are treated as weights, so if they do not add
    * up to 1 then they will be scaled appropriately. The `ws` need to be positive.
    * - It can have properties `f`, `min`, `max`, where `min` and `max` are numbers
    * and `f(i)` is a function. Then the distribution's values are meant to be the
    * sequential numbers (arithmetic sequence) from `min` to `max` with
    * an increment of 1, and `f(i)` gives the probability of the value `i`
    * (similar to the `ws`).
    *
    * `finite` returns an object representing the distribution, `dfinite` returns
    * the pmf, `pfinite` the cdf, `qfinite` the inverse cdf, and `rfinite` generates
    * random variates.
    *
    * Examples:
    * ```
    * finite({ min: 1, max: 10, f: function() { return 1; } })  // 10-sided fair die
    * finite({ xs: [0, 1], ws: [2, 1] }) // Bernoulli trial P(0) = 2/3
    * ```
    * @module distributions.finite
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var utils, rgen;

   utils = require('../utils');
   rgen = require('../rgen/rgen');

   /**
    * Returns an object representing a finite distribution, with properties
    * `d`, `p`, `q`, `r`.
    * ```
    * finite(o).d(x, logp)            // same as dfinite(o, logp)(x)
    * finite(o).p(x, lowerTail, logp) // same as pfinite(o, lowerTail, logp)(x)
    * finite(o).q(x, lowerTail, logp) // same as qfinite(o, lowerTail, logp)(x)
    * finite(o).r(n)                  // same as rfinite(o)(n)
    * ```
    * @memberof finite
   */
   function finite(o) {
      var sum, xs, ws, cumsLeft, cumsRight, i, finObj;

      o = typeof o.f === 'function' ? popObjFromFunc(o)
                                    : {
                                       xs: o.xs.slice(0),
                                       ws: o.ws.slice(0)
                                    };

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

   function popObjFromFunc(o) {
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
       * Returns the probability at `x` for the finite distribution
       * represented by object `o`.
       *
       * `logp` defaults to `false`; if `logp` is `true`, returns the
       * logarithm of the result.
       *
       * @fullName dfinite(o, logp)(x)
       * @memberof finite
       */
       function dfinite(o, logp) {
         var distr;

         distr = finite(o);
         return function(x) { return distr.d(x, logp); };
      },
      pfinite:
      /**
       * Evaluates the cumulative distribution function at `x`
       * for the finite distribution represented by object `o`:
       * $$\textrm{pfinite}(o)(x) = \sum_{k \leq x} dfinite(k)$$
       *
       * `lowerTail` defaults to `true`; if `lowerTail` is `false`, returns
       * the upper tail probability instead:
       * $$\textrm{pfinite}(o)(x) = \sum_{k > x} dfinite(k)$$
       *
       * `logp` defaults to `false`; if `logp` is `true`, returns the logarithm
       * of the result.
       *
       * @fullName pfinite(o, lowerTail, logp)(x)
       * @memberof finite
       */
       function pfinite(o, lowerTail, logp) {
         var distr;

         distr = finite(o);
         return function(q) { return distr.p(q, lowerTail, logp); };
      },
      qfinite:
      /**
       * Evaluates the quantile function for the finite distribution
       * specified by object `o`.
       * In general, for a discrete probability
       * distribution, the *quantile* is defined as the smallest domain value
       * `x` such that $F(x) \geq p$, where $F$ is the cumulative
       * distribution function.
       *
       * `p` is the desired probability ($0 \leq p \leq 1$).
       *
       * `lowerTail` defaults to `true`; if `lowerTail` is `false`, `p` is
       * interpreted as an upper tail probability.
       *
       * `logp` defaults to `false`; if `logp` is `true`, interprets `p` as
       * the logarithm of the desired probability.
       *
       * `qfinite` tries to invert `pfinite` but cannot be an exact inverse.
       * In particular, for `lowerTail = true`:
       *  - if asked for the smallest quantile for which the left area (<=) is 0, `qfinite` returns
       * `min`
       *  - if asked for the smallest quantile for which the left area (<=) is 1,
       * `qfinite` returns `max`.
       *
       * The edge cases for `lowerTail = false` are symmetrical to the preceding:
       * `qfinite` returns `min` or `max` for a right-tail area (>) of 1 or 0, respectively.
       *
       * @fullName qfinite(o, lowerTail, logp)(p)
       * @memberof finite
       */
       function qfinite(o, lowerTail, logp) {
         var distr;

         distr = finite(o);
         return function(p) { return distr.q(p, lowerTail, logp); };
      },
      rfinite:
      /**
       * Returns a random variate from the finite distribution
       * specified by object `o`.
       *
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
