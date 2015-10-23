(function(define) {'use strict';
define(function(require) {

   var C, twopi, sqrt2pi, Rational, pnorm, qnorm, pWrap;

   C = require('../constants');
   twopi = C.twopi;
   sqrt2pi = C.sqrt2pi;
   Rational = require('../rational');
   pWrap = require('../utils').pWrap;

   // density
   function dnorm(mu, sigma, logp) {
      var c;
      logp = logp === true;
      c = Math.log(twopi);
      return function(x) {
         var z, p;
         z = (x - mu) / sigma;
         p = -0.5 * (c + z * z);
         return logp ? p : Math.exp(p);
      };
   }

   // cdf
   pnorm = (function() {
      // segmented rational approximations
      var small, medium, large;

      function addExpTerm(x, rational) {
         return Math.exp(-0.5 * x * x) * rational;
      }

      // Computations taken from GSL (gauss.c)
      small = Rational.new([
         0.065682337918207449113,
         2.2352520354606839287,
         161.02823106855587881,
         1067.6894854603709582,
         18154.981253343561249
      ], [
         1,
         47.20258190468824187,
         976.09855173777669322,
         10260.932208618978205,
         45507.789335026729956
      ]);

      medium = Rational.new([
         1.0765576773720192317e-8,
         0.39894151208813466764,
         8.8831497943883759412,
         93.506656132177855979,
         597.27027639480026226,
         2494.5375852903726711,
         6848.1904505362823326,
         11602.651437647350124,
         9842.7148383839780218
      ], [
         1,
         22.266688044328115691,
         235.38790178262499861,
         1519.377599407554805,
         6485.558298266760755,
         18615.571640885098091,
         34900.952721145977266,
         38912.003286093271411,
         19685.429676859990727
      ]);

      large = Rational.new([
         0.02307344176494017303,
         0.21589853405795699,
         0.1274011611602473639,
         0.022235277870649807,
         0.001421619193227893466,
         2.9112874951168792e-5
      ], [
         1,
         1.28426009614491121,
         0.468238212480865118,
         0.0659881378689285515,
         0.00378239633202758244,
         7.29751555083966205e-5
      ]);

      return function pnorm(mu, sigma, lowerTail, logp) {
         lowerTail = lowerTail !== false;
         logp = logp === true;
         /* eslint-disable complexity */
         return function(x) {
            var z, absz, oneOverzsq, r, ret, lower;

            lower = lowerTail;
            z = (x - mu) / sigma;
            if (z > 5 && !lower) { /* danger zone for upper */
               z = -z;
               lower = true;
            }
            absz = Math.abs(z);

            if (absz < 0.66291) {
               ret = 0.5 + z * small.evalAt(z * z);
            } else if (absz < 5.656854) {
               r = addExpTerm(z, medium.evalAt(absz));
               ret = z > 0 ? 1 - r : r;
            } else if (z > 8.572) {
               ret = 1;
            } else if (z < -37.519) {
               ret = 0;
            } else {
               oneOverzsq = 1 / (z * z);
               r = oneOverzsq * large.evalAt(oneOverzsq);
               r = addExpTerm(z, (1 / sqrt2pi - r) / absz);
               ret = z > 0 ? 1 - r : r;
            }
            if (!lower) { ret = 1 - ret; }
            if (logp) { ret = Math.log(ret); }
            return ret;
         };
         /* eslint-enable complexity */
      };
   }());

   // inverse cdf
   qnorm = (function() {
   // segmented rational approximations
      var small, medium, large;

      small = Rational.new([
         2509.0809287301226727,
         33430.575583588128105,
         67265.770927008700853,
         45921.953931549871457,
         13731.693765509461125,
         1971.5909503065514427,
         133.14166789178437745,
         3.387132872796366608
      ], [
         5226.495278852854561,
         28729.085735721942674,
         39307.89580009271061,
         21213.794301586595867,
         5394.1960214247511077,
         687.1870074920579083,
         42.313330701600911252,
         1.0
      ]);

      medium = Rational.new([
         7.7454501427834140764e-4,
         0.0227238449892691845833,
         0.24178072517745061177,
         1.27045825245236838258,
         3.64784832476320460504,
         5.7694972214606914055,
         4.6303378461565452959,
         1.42343711074968357734
      ], [
         1.05075007164441684324e-9,
         5.475938084995344946e-4,
         0.0151986665636164571966,
         0.14810397642748007459,
         0.68976733498510000455,
         1.6763848301838038494,
         2.05319162663775882187,
         1.0
      ]);

      large = Rational.new([
         2.01033439929228813265e-7,
         2.71155556874348757815e-5,
         0.0012426609473880784386,
         0.026532189526576123093,
         0.29656057182850489123,
         1.7848265399172913358,
         5.4637849111641143699,
         6.6579046435011037772
      ], [
         2.04426310338993978564e-15,
         1.4215117583164458887e-7,
         1.8463183175100546818e-5,
         7.868691311456132591e-4,
         0.0148753612908506148525,
         0.13692988092273580531,
         0.59983220655588793769,
         1.0
      ]);

      return function qnorm(mu, sigma, lowerTail, logp) {
         lowerTail = lowerTail !== false;
         logp = logp === true;
         return pWrap(lowerTail, logp, function(p) {
            var dp, minp, r;

            dp = p - 0.5;
            minp = p < 0.5 ? p : 1 - p;

            if (p === 1) { return Infinity; }
            if (p === 0) { return -Infinity; }
            if (Math.abs(dp) <= 0.425) {
               r = dp * small.evalAt(0.180625 - dp * dp);
            } else {
               r = Math.sqrt(-Math.log(minp));
               r = r <= 5 ? medium.evalAt(r - 1.6) : large.evalAt(r - 5);
               r = p < 0.5 ? -r : r;
            }
            return mu + sigma * r;
         });
      };
   }());

   return {
      dnorm: dnorm,
      pnorm: pnorm,
      qnorm: qnorm
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
