(function(define) {
'use strict';
define(function(require) {

   var utils, log1p, expm1;

   log1p = require('./basicFunc/log1p').log1p;
   expm1 = require('./basicFunc/expm1').expm1;

   utils = {
      /* mixin */
      mixin: function mixin(target) {
         Array.prototype.slice.call(arguments, 1)
            .forEach(function(source) {
               Object.keys(source).forEach(function(key) {
                  target[key] = source[key];
               });
            });
         return target;
      },
      /* Maximum number of steps for `repeat` iterations. */
      maxSteps: 1000,
      /* repeat, for carrying out repeated improvement until stopping
       * condition is satisfied, or until `utils.maxSteps` have been
       * performed.
       *
       * `init` is the initial value,
       * `step` is a function which returns the new value, and
       * `stop` is a non-negative integer (number of steps to take),
       *  a function of no arguments returning a boolean, or missing.
       *  If `stop` is omitted, the default behavior is using
       * `utils.relativelyCloseTo`.
       *  Iteration will stop after `utils.maxSteps`, regardless.
       */
      repeat: function(init, step, stop) {
         var prev, curr, done, reps;

         curr = init;
         reps = utils.maxSteps;
         if (stop === 0) { return curr; }
         done = typeof stop === 'function' ?
                stop
              : stop > 0 ? function() { stop -= 1; return stop === 0; }
                         : function() { return utils.relativelyCloseTo(curr, prev, 1e-13); };
         while (reps > 0 && !done() && !isNaN(curr)) {
            reps -= 1;
            prev = curr;
            curr = step();
         }
         if (isNaN(curr)) { console.log('repeat returned NaN'); }
         if (reps === 0) {
            console.log('repeat terminated after too many repetitions');
         }
         return curr;
      },
      /* series
       * - `f` is a function two arguments, `index` and `prev` (previous
       *     value).  It will be called starting at `index` 0.
       * - `stop` is the same as in `repeat`.
       * Add up the terms generated by `f` for index >= 0; stopping conditions
       * are the same as for `repeat`.
       */
      series: function(f, stop) {
         var sum, curr, i;

         i = 0;
         curr = f(0);
         sum = curr;
         return utils.repeat(sum, function() {
            i += 1;
            curr = f(i, curr);
            sum += curr;
            return sum;
         }, stop);
      },
      /* contFrac
       * A continued fraction has the form
       * a0 + (b1 / (a1 + (b2 / (a2 + ...)))
       * where a is indexed from 0 and b is indexed from 1.    *
       * - `a` and `b` are functions with parameters `i` and `v` for generating
       * the ai and bi sequences; each time, the function returns the ith term
       * of the sequence based on the previous value, `v`.
       * - `stop` is the same as in `repeat`.
       *
       * Return a value for the indicated continued fraction, to a precision
       * as determined by `stop`.
       */
      contFrac: function(a, b, stop) {
         var A, A1, A2, // A2 is "two previous" to A, A1 is "one previous"
             B, B1, B2, an, bn, i;

         an = a(0);
         A = an;
         A1 = 1;
         bn = 1;
         B = bn;
         B1 = 0;
         i = 0;
         return utils.repeat(A / B, function() {
            A2 = A1;
            A1 = A;
            B2 = B1;
            B1 = B;
            i += 1;
            an = a(i, an);
            bn = b(i, bn);
            A = A1 * an + A2 * bn;
            B = B1 * an + B2 * bn;
            return A / B;
         }, stop);
      },
      /* Binary search for solution to f(x) = y (provided f is strictly monotone)
       *
       * a, b are optional bounds. If not provided,
       * the function will search for them using x = +/- 2^n
       *
       * Search will terminate when:
       * - Provided a, b values do not satisfy f(a) <= y <= f(b) (error)
       * - The search for a, b produces unreasonably large endpoints (error)
       * - The relative size of the search interval becomes sufficiently small
       *      relative to its midpoint.
       * - utils.maxSteps is reached (error)
       * - utils.relativelyCloseTo (f(x), y) < utils.precision
       */
       /* eslint-disable complexity */
      binSearchSolve: function binSearchSolve(f, y, a, b) {
         var mid, vmid, niters;

         a = inferEndpoint(f, y, true, a);
         b = inferEndpoint(f, y, false, b);

         if (f(a) > f(b)) {
            return binSearchSolve(function(x) { return -f(x); }, -y, a, b);
         }

         if (utils.relativelyCloseTo(f(a), y)) { return a; }
         if (utils.relativelyCloseTo(f(b), y)) { return b; }

         if (f(a) > y || f(b) < y) {
            throw new Error('Binary search: Desired value not between endpoint values.');
         }
         niters = 0;
         while (niters < utils.maxSteps) {
            mid = (a + b) / 2;
            niters += 1;
            vmid = f(mid);
            if (utils.relativelyCloseTo(a, b) &&
                utils.relativelyCloseTo(vmid, y)) {
              return mid;
            }
            if (vmid < y) {
               a = mid;
            } else {
               b = mid;
            }
         }
         throw new Error('Binary search: Too many iterations');
      }, /* eslint-enable complexity */
      /* eslint-disable complexity */
      /*
       *  Helper for inverting a discrete CDF for integer values.
       *  Searches for a quantile x in the range [`min`, `max`] such
       *  that x is the smallest quantile for the provided CDF `f`
       *  with left-tail probability >= `p`.
       *  `x` is an initial estimate for the desired quantile.
       *  When `p` is 0, return min.
       *  When `p` is 1, return max.
       *  Precondition:  `0 <= p <= 1`.
       *  Precondition:  `f(x)` is 0 when `x < min` and 1 when `x > max`.
       */
      discInvCdf: function discInvCdf(min, max, x, p, f) {
         var incr;

         if (p === 0) { return min; }
         if (p === 1) { x = max; }
         if (x > max) { x = max; }
         if (x < min) { x = min; }

         incr = Math.min(x - min, max - x);
         incr = Math.max(Math.floor(0.001 * incr), 1);
         while (incr > 1) {
            if (f(x) < p) { // x is too small
               x += incr;
               if (f(x) >= p) {
                  incr = Math.floor(incr / 2);
               }
            } else {
               // x is probably too large
               x -= incr;
               if (f(x) < p) {
                  incr = Math.floor(incr / 2);
               }
            }
         }
         if (x < min) { x = min; }
         if (x > max) { x = max; }
         // incr is 1
         /* eslint-disable no-constant-condition */
         while (true) {
            if (x === min) {
               if (f(x) >= p || utils.relativelyCloseTo(f(x), p, 1e-14)) { return x; }
               x += 1;  // go right
            } else if (f(x - 1) >= p ||
               utils.relativelyCloseTo(f(x - 1), p, 1e-14)) {
               if (Math.abs(f(x) - p) < Math.abs(f(x - 1) - p)) {
                  return x;  // don't go left, x is better!
               }
               x -= 1; // go left
            } else {
               // return x, or go right?
               if (x === max || f(x) >= p || utils.relativelyCloseTo(f(x), p, 1e-14)) {
                  return x;
               }
               x += 1;  // go right
            }
         }
         /* eslint-enable no-constant-condition */
      },
       /* eslint-enable complexity */

      /* precision used by relativelyCloseTo */
      precision: 1e-10,
      /* relativelyCloseTo returns a boolean indicating whether x, x0 are
       * relatively close to each other as specified by the precision `delta`.
       * If `delta` is not provided, `utils.precision` is used instead.
       *
       * Special cases for NaN, Infinity, -Infinity
       */
      relativelyCloseTo: function(x, x0, delta) {
         var absMax;

         delta = delta || utils.precision;
         absMax = Math.max(Math.abs(x0), Math.abs(x));
         if (isNaN(absMax)) { return isNaN(x) && isNaN(x0); } /* both NaN */
         if (absMax === Infinity) { return x0 === x; }
         if (absMax === 0) { return true; }
         if (x === 0) { return utils.isEssentiallyZero(x0); }
         if (x0 === 0) { return utils.isEssentiallyZero(x); }
         return Math.abs(x - x0) / absMax < delta;
      },
      /* precision used by relativelyCloseTo */
      /* isEssentiallyZero returns a boolean indicating whether x is very
       * close to zero (less than e-100 in absolute value).
       */
      isEssentiallyZero: function(x) {
         return Math.abs(x) < 1e-100;
      },
      /*
       * If lx = log(x), ly = log(y), calculates log(x + y)
       */
      logspaceAdd: function(lx, ly) {
         return Math.max(lx, ly) + log1p(Math.exp(-Math.abs(lx - ly)));
      },
      /*
       * Given a lower tail probability value, it "adjusts" it to
       * requested lowerTail and logp settings.
       *
       * It expects p in the range [0, 1]
       */
      adjustLower: function(p, lowerTail, logp) {
        if (p === 1) {
          p = 0;
          lowerTail = !lowerTail;
        }
        if (p === 0) {
          return lowerTail ? logp ? -Infinity : 0
                           : logp ? 0 : 1;
        }
        return lowerTail ? logp ? Math.log(p) : p
                         : logp ? log1p(-p) : 1 - p;
      },
      /*
       * Given a upper tail probability value, it "adjusts" it to
       * requested lowerTail and logp settings.
       *
       * It expects q in the range (0, 1).
       * Use `adjustLower` for q = 0 or 1.
       */
      adjustUpper: function(q, lowerTail, logp) {
        return lowerTail ? logp ? log1p(-q) : 1 - q
                         : logp ? Math.log(q) : q;
      },
      /*
       * Given a p-value in lowerTail/logp space, returns an object
       * of the logged lowerTail/upperTail probabilities p, q
       */
      logProbs: function logProbs(prob, lowerTail, logp) {
         return logp ? lowerTail ? { p: prob, q: Math.log(-expm1(prob)) }
                                 : { p: Math.log(-expm1(prob)), q: prob }
                     : lowerTail ? { p: Math.log(prob), q: log1p(-prob) }
                                 : { p: log1p(-prob), q: Math.log(prob) };
      },
      /*
       * Given a p-value in lowerTail/logp space, returns an object
       * of the unlogged lowerTail/upperTail probabilities p, q
       */
      trueProbs: function trueProbs(prob, lowerTail, logp) {
         return logp ? lowerTail ? { p: Math.exp(prob), q: -expm1(prob) }
                                 : { p: -expm1(prob), q: Math.exp(prob) }
                     : lowerTail ? { p: prob, q: 1 - prob }
                                 : { p: 1 - prob, q: prob };
      },
      /*
       * Given a function 'f(prob)', returns a function of the probability
       * which has been adjusted in the space specified by 'lowerTail' and
       * 'logp'.  The new function wraps f, ensuring that f is called on
       * an actual left-tail probability (unlogged, complemented if
       * necessary).
       */
      pWrap: function pWrap(lowerTail, logp, f) {
         return function(prob) {
            var probs;

            probs = utils.trueProbs(prob, lowerTail, logp);

            return probs.p >= 0 && probs.q >= 0 ? f(probs) : NaN;
         };
      },

      /*
       * qhelper
       * helper for quantile functions, to handle out-of-range cases
       *
       * Takes `logp` and a function `f`.
       * Returns a function of `p` which checks that `p` is a valid probability
       * in the logp space before passing `p` to `f`.
       *
       * The returned function returns NaN if the probability represented by `p` is out
       * of range.
       */
      qhelper: function qhelper(logp, f) {
         return logp ? function(p) { return p <= 0 ? f(p) : NaN; }
                     : function(p) { return p >= 0 && p <= 1 ? f(p) : NaN; };
      }
   };

   // lower = true means looking for lower endpoint
   //       = false means looking for upper endpoint
   // if "e" is provided it is just returned
   function inferEndpoint(f, y, lower, e) {
      var comp; // comparison of f(e) and f(2*e)

      if (e != null) { return e; }

      e = lower ? -1 : 1;
      comp = f(e) < f(2 * e);

      while (f(e) > y !== comp) {
         e *= 2;
         if (-e > 1e200) {
            throw new Error('Binary search: Cannot find ' + (lower ? 'lower' : 'upper') + ' endpoint.');
         }
      }
      return e;
   }

   return utils;

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
