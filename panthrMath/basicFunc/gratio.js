(function(define) {
'use strict';
/* eslint-disable max-statements */
define(function(require) {

   /*
    * Implementation of the incomplete gamma function and its
    * complement.
    * See: "Computation of the Incomplete Gamma Function Ratios and their
    * Inverse", DiDonato and Morris
    *
    * P(a, x), Q(a, x) are the ratios, where P + Q = 1.
    * Domain: a >= 0, x >= 0, a + x != 0.
    */

   var phi, gam1, repeat, series, contFrac, cf, BIG, x0, e0, t, logroot, r,
       Rational, Polynomial, gamma, lgamma, erf, erfc, expm1,
       smallP, smallQ, mediumP, mediumQ, bigP, bigQ,
       eulerGamma, funS;

   Rational = require('../rational');
   Polynomial = require('../polynomial');
   gamma = require('./lgamma').gamma;
   lgamma = require('./lgamma').lgamma;
   erf = require('./erf').erf;
   erfc = require('./erf').erfc;
   repeat = require('../utils').repeat;
   series = require('../utils').series;
   contFrac = require('../utils').contFrac;
   expm1 = require('./expm1').expm1;
   eulerGamma = require('../constants').eulerGamma;
   gam1 = require('./gam1').gam1;
   phi = require('./phi').phi;

   logroot = Math.log(Math.sqrt(0.765));

   // See table 7
   BIG = 20;
   x0 = 31;
   e0 = 0.25e-3;

   // Function T(a, lambda) (Formula 18)
   t = (function() {
      var cs;

      cs = [ // D0 - degree 13
         new Polynomial([
            -0.438203601845335318655297462245e-8,
             0.102618097842403080425739573227e-7,
             0.670785354340149858036939710030e-8,
            -0.176659527368260793043600542457e-6,
             0.829671134095308600501624213166e-6,
            -0.185406221071515996070179883623e-5,
            -0.218544851067999216147364295512e-5,
             0.391926317852243778169704095630e-4,
            -0.178755144032921810699588477366e-3,
             0.352733686067019400352733686067e-3,
             0.115740740740740740740740740741e-2,
            -0.148148148148148148148148148148e-1,
             0.833333333333333333333333333333e-1,
            -0.333333333333333333333333333333
         ]), // D1 - degree 12
         new Polynomial([
             0.119516285997781473243076536700e-7,
            -0.575254560351770496402194531835e-7,
             0.137863344691572095931187533077e-6,
             0.464712780280743434226135033939e-8,
            -0.161209008945634460037752218822e-5,
             0.764916091608111008463742149809e-5,
            -0.180985503344899778370285914868e-4,
            -0.401877572016460905349794238683e-6,
             0.205761316872427983539094650206e-3,
            -0.990226337448559670781893004115e-3,
             0.264550264550264550264550264550e-2,
            -0.347222222222222222222222222222e-2,
            -0.185185185185185185185185185185e-2
         ]), // D2 - degree 10
         new Polynomial([
             0.142806142060642417915846008823e-6,
            -0.629899213838005502290672234278e-6,
             0.137219573090629332055943852926e-5,
             0.342357873409613807419020039047e-7,
            -0.127606351886187277133779191392e-4,
             0.529234488291201254164217127180e-4,
            -0.107366532263651605215391223622e-3,
             0.200938786008230452674897119342e-5,
             0.771604938271604938271604938272e-3,
            -0.268132716049382716049382716049e-2,
             0.413359788359788359788359788360e-2
         ]), // D3 - degree 8
         new Polynomial([
             0.142309007324358839145518944706e-5,
            -0.567495282699159656749963105702e-5,
             0.110826541153473023614770299727e-4,
            -0.239650511386729665193314027333e-6,
            -0.756180167188397641072538191880e-4,
             0.267720632062838852962309752433e-3,
            -0.469189494395255712128140111679e-3,
             0.229472093621399176954732510288e-3,
             0.649434156378600823045267489712e-3
         ]), // D4 - degree 6
         new Polynomial([
             0.113757269706784190980552042886e-4,
            -0.396836504717943466443123507595e-4,
             0.664149821546512218665853782452e-4,
            -0.146384525788434181781232535691e-5,
            -0.299072480303190179733389609933e-3,
             0.784039221720066627474034881442e-3,
            -0.861888290916711698604702719929e-3
         ]), // D5 - degree 4
         new Polynomial([
             0.679778047793720783881640176604e-4,
            -0.199325705161888477003360405281e-3,
             0.277275324495939207873364251965e-3,
            -0.697281375836585777429398828576e-4,
            -0.336798553366358150308767592718e-3
         ]), // D6 - degree 2
         new Polynomial([
             0.270878209671804482771279183488e-3,
            -0.592166437353693882864836225604e-3,
             0.531307936463992223165748542978e-3
         ])
      ];

      return function(a, lambda) {
         var z, atomk; // atomk: a^{-k}

         atomk = a;
         z = (lambda >= 1 ? 1 : -1) * Math.sqrt(2 * phi(lambda));

         return cs.reduce(function(sum, poly) {
            atomk /= a;
            return sum + poly.evalAt(z) * atomk;
         }, 0);
      };
   }());

   function e(y) {
      return 0.5 - (1 - y / 3) * Math.sqrt(y / Math.PI);
   }

   function j(a, x) {
      var v;

      v = -x;
      return -a * series(function(i) {
         if (i === 0) { return 0; }
         if (i === 1) { return -x / (a + 1); }
         v = v * -x / i;
         return v / (a + i);
      });
   }

   r = (function() {
      // From equation 3
      function delta(a) {
         return lgamma(a) - (a - 0.5) * Math.log(a) + a - 0.5 * Math.log(2 * Math.PI);
      }

      return function(a, x) {
         return a <= 20 ? Math.exp(-x) * Math.pow(x, a) / gamma(a)
                        : Math.sqrt(a / (2 * Math.PI)) *
                          Math.exp(-a * phi(x / a) - delta(a));
      };
   }());

   // cf(a, x) computes the continued fraction in equation 11
   cf = function(a, x) {
      function as(i, v) {
         return i === 0 ? 0
                        : i % 2 === 1 ? x
                                      : 1;
      }
      function bs(i, v) {
         if (i % 2 === 1) {
            return i === 1 ? 1 : (i - 1) / 2;
         }
         return i / 2 - a;
      }
      return contFrac(as, bs);
   };

   smallP = function(a) {
      var gamm;

      gamm = gamma(a + 1);
      return function(x) {
         var alpha;

         if (x < 1.1) {
            alpha = x < 0.5 ? logroot / Math.log(x)
                            : x / 2.59;
            return a >= alpha ? Math.pow(x, a) * (1 - j(a, x)) / gamm // Formula 9
                              : 1 - smallQ(a)(x);
         }
         return 1 - smallQ(a)(x);
      };
   };

   smallQ = function(a) {
      var gamm, ha;

      gamm = gamma(a + 1);
      ha = gam1(a);
      return function(x) {
         var alpha;

         if (x < 1.1) {
            alpha = x < 0.5 ? logroot / Math.log(x)
                            : x / 2.59;
            return a >= alpha ? 1 - smallP(a)(x)   // Formula 10
                              : (Math.pow(x, a) * j(a, x) - expm1(a * Math.log(x))
                                ) / gamm - ha;
         } // Formula 11
         return r(a, x) * cf(a, x);
      };
   };

   mediumP = function(a) {
      return function(x) {
         if ((a > x || x >= x0 || a !== Math.round(2 * a)) &&
             x <= Math.max(a, Math.LN10)) {
            // Formula 15, for P
            return r(a, x) / a * series(function(n, v) {
               return n === 0 ? 1 : v * x / (a + n);
            });
         }
         return 1 - mediumQ(a)(x);
      };
   };

   mediumQ = function(a) {
      return function(x) {
         if (a <= x && x < x0 && 2 * a === Math.round(2 * a)) {
            // Formulas 14
            if (a === Math.floor(a)) { // a integer
               return Math.exp(-x) * series(function(i, v) {
                  return i === 0 ? 1 : v * x / i;
               }, a);
            } // a half-integer
            return erfc(Math.sqrt(x)) +
                   Math.exp(-x) / Math.sqrt(Math.PI * x) *
                     series(function(i, v) {
                        return i === 0 ? 0
                                       : i === 1 ? x / 0.5
                                                 : v * x / (i - 0.5);

                     }, Math.round(a + 0.5));
         }
         if (x <= Math.max(a, Math.LN10)) {
            // Formula 15, for P
            return 1 - mediumP(a)(x);
         }
         if (x < x0) { // Formula 11, for Q
            return r(a, x) * cf(a, x);
         }
         // Formula 16, for Q
         return r(a, x) / x * series(function(n, v) {
            return n === 0 ? 1 : v * (a - n) / x;
         });
      };
   };

   bigP = function(a) {
      return function(x) {
         var lambda, sigma, y;

         lambda = x / a;
         sigma = Math.abs(1 - lambda);
         y = a * phi(lambda);

         if (sigma <= e0 / Math.sqrt(a)) {
            // Formula 19, for P
            if (lambda > 1) { return 1 - bigQ(a)(x); }
            return e(y) - (1 - y) / Math.sqrt(2 * Math.PI * a) * t(a, lambda);
         }
         if (sigma <= 0.4) {
            // Formula 17, for P
            if (lambda > 1) { return 1 - bigQ(a)(x); }
            return 0.5 * erfc(Math.sqrt(y)) -
               Math.exp(-y) / Math.sqrt(2 * Math.PI * a) * t(a, lambda);
         }
         if (x <= Math.max(a, Math.LN10)) {
            // Formula 15, for P
            return r(a, x) / a * series(function(n, v) {
               return n === 0 ? 1 : v * x / (a + n);
            });
         }
         return 1 - bigQ(a)(x);
      };
   };

   bigQ = function(a) {
      return function(x) {
         var lambda, sigma, y;

         lambda = x / a;
         sigma = Math.abs(1 - lambda);
         y = a * phi(lambda);

         if (sigma <= e0 / Math.sqrt(a)) {
            // Formula 19, for Q
            if (lambda <= 1) { return 1 - bigP(a)(x); }
            return e(y) + (1 - y) / Math.sqrt(2 * Math.PI * a) * t(a, lambda);
         }
         if (sigma <= 0.4) {
            // Formula 17, for Q
            if (lambda <= 1) {
               return 1 - bigP(a)(x);
            }
            return 0.5 * erfc(Math.sqrt(y)) +
               Math.exp(-y) / Math.sqrt(2 * Math.PI * a) * t(a, lambda);
         }
         if (x <= Math.max(a, Math.LN10)) {
            // Formula 15, for P
            return 1 - bigP(a)(x);
         }
         if (x < x0) { // Formula 11, for Q
            return r(a, x) * cf(a, x);
         }
         // Formula 16, for Q
         return r(a, x) / x * series(function(n, v) {
            return n === 0 ? 1 : v * (a - n) / x;
         });

      };
   };

   /**
    * Returns the incomplete gamma ratio:
    * $$P(a,x) = \frac{1}{\Gamma(a)} \int_0^x e^{-t} t^{a-1} dt$$
    *
    * `a` is the shape parameter and $0 < x < \infty$.
    *
    * See also: `pgamma`
    *
    * Based on:  *Computation of the Incomplete Gamma Function Ratios
    * and their Inverse*, by DiDonato and Morris, 1986
    * @fullName gratio(a)(x)
    * @memberof basicFunc
    */
   function gratio(a) {
      if (a <= 0) { return function(x) { return NaN; }; }
      if (a === 0.5) { return function(x) { return erf(Math.sqrt(x)); }; }
      if (a < 1) { return smallP(a); }
      if (a < BIG) { return mediumP(a); }
      return bigP(a);
   }

   /**
    * Returns the complement of the incomplete gamma ratio:
    * $$Q(a,x) = \frac{1}{\Gamma(a)} \int_x^\infty e^{-t} t^{a-1} dt$$
    *
    * `a` is the shape parameter and $0 < x < \infty$.
    *
    * See also: `pgamma`
    *
    * Based on:  *Computation of the Incomplete Gamma Function Ratios
    * and their Inverse*, by DiDonato and Morris, 1986
    * @fullName gratioc(a)(x)
    * @memberof basicFunc
    */
   function gratioc(a) {
      if (a <= 0) { return function(x) { return NaN; }; }
      if (a === 0.5) { return function(x) { return erfc(Math.sqrt(x)); }; }
      if (a < 1) { return smallQ(a); }
      if (a < BIG) { return mediumQ(a); }
      return bigQ(a);
   }

   // INVERSE GAMMA

   // funS implements formula 32
   // a function of p and q for calculating s
   funS = (function() {
      var ratFun;

      ratFun = Rational.new([
         0.213623493715853, 4.28342155967104,
         11.6616720288968, 3.31125922108741
      ], [
         0.0361170810188420, 1.27364489782223,
         6.40691597760039, 6.61053765625462, 1
      ]);
      return function(p, q) {
         var temp;

         temp = Math.sqrt(-2 * Math.log(p < 0.5 ? p : q));
         return (p < 0.5 ? -1 : 1) * (temp - ratFun.evalAt(temp));
      };
   }());

   // Use formulas 21 through 25 from DiDonato to get initial
   // estimate for x when a < 1.
   function findx0SmallA(a, p, q) {
      var B, u, temp;

      B = q * gamma(a);
      if (B > 0.6 || B >= 0.45 && a >= 0.3) {
         // use 21
         u = B * q > 1e-8 ? Math.pow(p * gamma(a + 1), 1 / a)
                          : Math.exp(-q / a - eulerGamma);
         return u / (1 - u / (a + 1));
      }
      if (a < 0.3 && B >= 0.35 && B <= 0.6) {
         // use 22
         temp = Math.exp(-eulerGamma - B);
         u = temp * Math.exp(temp);
         return temp * Math.exp(u);
      }
      temp = -Math.log(B);
      u = temp - (1 - a) * Math.log(temp);
      if (B >= 0.15) {
         // use 23, with temp for y and u for v
         return temp - (1 - a) * Math.log(u) -
                  Math.log(1 + (1 - a) / (1 + u));
      }
      if (B > 0.01) {
         // use 24, with temp for y and u for v
         return temp - (1 - a) * Math.log(u) - Math.log(
            (u * u + 2 * (3 - a) * u + (2 - a) * (3 - a)) /
            (u * u + (5 - a) * u + 2));
      }
      // use 25, where c1 + ... + c5/y^4 is treated as a polynomial in y^-1
      // using u for c1, temp for y
      return findx0TinyB(a, temp);
   }

   // implements equation #25
   function findx0TinyB(a, y) {
      var c1;

      c1 = (a - 1) * Math.log(y);
      return y + Polynomial.new([
         // c5
         (a - 1) * Polynomial.new([
            -0.25, (11 * a - 17) / 6, -3 * a * a + 13 * a - 13,
            (2 * a * a * a - 25 * a * a + 72 * a - 61) / 2,
            (25 * a * a * a - 195 * a * a + 477 * a - 379) / 12
         ]).evalAt(c1),
         // c4
         (a - 1) * Polynomial.new([
            1 / 3, (3 * a - 5) / 2, a * a - 6 * a + 7,
            (11 * a * a - 46 * a + 47) / 6
         ]).evalAt(c1),
         // c3
         (a - 1) * (-0.5 * c1 * c1 + (a - 2) * c1 + (3 * a - 5) / 2),
         // c2
         (a - 1) * (1 + c1),
         // c1
         c1
      ]).evalAt(1 / y);
   }

   // Use formulas 25 and 31 through 36 from DiDonato to get initial
   // estimate for x when a > 1.
   /* eslint-disable complexity */
   function findx0BigA(a, p, q) {
      var w, s, D, B, u, z, zbar, logpg;

      // snTerm -- see formula 34
      logpg = Math.log(p * gamma(a + 1));
      function snTerm(x) {
         return function(n, v) {
            return n === 0 ? 1 : v * x / (a + n);
         };
      }
      // See formula 34, where `n` is off by 1
      function f(x, n) {
         return Math.exp((logpg + x - Math.log(series(snTerm(x), n + 1))) / a);
      }
      B = q * gamma(a);
      D = Math.max(2, a * (a - 1));
      s = funS(p, q);
      w = a + s * Math.sqrt(a) + (s * s - 1) / 3 +
         (s * s * s - 7 * s) / (36 * Math.sqrt(a)) -
         (3 * s * s * s * s + 7 * s * s - 16) / (810 * a) +
         (9 * s * s * s * s * s + 256 * s * s * s - 433 * s) /
            (38880 * a * Math.sqrt(a));
      if (a >= 500 && Math.abs(1 - w / a) < 1e-6) { return w; }
      if (p > 0.5) {
         if (w < 3 * a) { return w; }
         if (B < Math.pow(10, -D)) { return findx0TinyB(a, -Math.log(B)); }
         u = -Math.log(B) + (a - 1) * Math.log(w) -
               Math.log(1 + (1 - a) / (1 + w));
         // formula 33
         return -Math.log(B) + (a - 1) * Math.log(u) -
               Math.log(1 + (1 - a) / (1 + u));
      }
      // handle p <= 0.5
      z = w > 0.15 * (a + 1) ? w : f(f(f(f(w, 1), 3), 3), 4);
      if (z <= 0.01 * (a + 1) || z > 0.7 * (a + 1)) { return z; }
      /* eslint-disable no-extra-parens */
      return (function() {
         // find N for formula 36
         var N, term;

         N = 0;
         term = 1;
         do {
            N += 1;
            term = snTerm(z)(N, term);
         } while (term > 1e-4);
         zbar = f(z, N + 1);
         return zbar * (1 - (a * Math.log(zbar) - zbar - logpg +
                  Math.log(series(snTerm(z), N + 1))) / (a - zbar));
      }());
      /* eslint-enable no-extra-parens */
   }
   /* eslint-enable complexity */

   // Use Schroeder or Newton-Raphson to do one step of the iteration,
   // formulas 37 and 38.
   function step(x, a, p, q) {
      var temp, w;

      temp = p <= 0.5 ? gratio(a)(x) - p : q - gratioc(a)(x);
      temp /= r(a, x);
      w = (a - 1 - x) / 2;
      if (Math.max(Math.abs(temp), Math.abs(w * temp)) <= 0.1) {
         return x * (1 - (temp + w * temp * temp));
      }
      return x * (1 - temp);
   }

    /**
     * Returns the inverse function
     * of the incomplete gamma ratio (see: `gratio`).
     *
     * `a` is the shape parameter and `p` is the desired tail probability.
     * `lower` is a boolean which determines whether `p` is a left-tail
     * probability or not (defaults to `true`).
     *
     * See also: `qgamma`
     *
     * Based on:  *Computation of the Incomplete Gamma Function Ratios
     * and their Inverse*, by DiDonato and Morris, 1986
     * @fullName gaminv(a)(p, lower)
     * @memberof basicFunc
     */
   function gaminv(a) {
      if (a <= 0) { return function(p) { return NaN; }; }
      /* eslint-disable complexity */
      return function(p, lower) {
         var q, x;

         if (lower === false) {
            q = p;
            p = 1 - q;
         } else {
            q = 1 - p;
         }
         if (p === 0) { return 0; }
         if (p === 1) { return Infinity; }
         if (p < 0 || p > 1) { return NaN; }

         if (a === 1) { return -Math.log(q); }
         x = a < 1 ? findx0SmallA(a, p, q) : findx0BigA(a, p, q);
         return repeat(x, function() {
            x = step(x, a, p, q);
            return x;
         });
      };
      /* eslint-enable complexity */
   }

   /* eslint-disable complexity */

  /**
   * Computes the scaled complement of the incomplete gamma ratio function
   * (see: `gratio`): $$\textrm{grat}_r(a,x,r) =  \frac{Q(a,x)}{r}$$
   * where
   *               $Q(a,x)$ is the upper-tail incomplete gamma function
   * (see: `gratioc`) and  $r = \frac{e^{-x} x^a}{\Gamma(a)}$
   *
   * Meant to be used internally.
   *
   * Based on:   *Computation of the Incomplete Gamma Function Ratios
   * and their Inverse*, by DiDonato and Morris, 1986
   * @memberof basicFunc
   */
   function gratR(a, x, logr) {
      var ser, z, l, h, rinv;

      rinv = Math.exp(-logr);
      /* R's TOMS code handles this slightly differently (1214 in toms708.c) */
      if (x === 0) { return rinv; }
      if (a === 0) { return 0; }
      if (a === 0.5) {
         return x < 0.25 ? (1 - erf(Math.sqrt(x))) * rinv
                         : Math.exp(x) * erfc(Math.sqrt(x)) / Math.sqrt(x) * Math.sqrt(Math.PI);
      }
      z = a * Math.log(x);
      l = expm1(z);
      h = gam1(a);
      if (x < 1.1) {
         ser = a * x * series((function(c, v) {
               return function(i) {
                  c += 1;
                  if (i > 0) { v = v * -x / (i + 1); }
                  return v / c;
               };
            }(a, 1))
         );
         if (x >= 0.25 && a < x / 2.59 || z > -0.13394) {
            return Math.max(0, ((Math.exp(z) * ser - l) * (1 + h) - h) * rinv);
         }
         return (1 - Math.exp(z) * (1 + h) * (1 - ser)) * rinv;
      }
      return cf(a, x);
   }

   return {
      gratio: gratio,
      gratioc: gratioc,
      gaminv: gaminv,
      gratR: gratR
   };
   /* eslint-enable complexity */

});
/* eslint-enable max-statements */
}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
