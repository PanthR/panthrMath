(function(define) {
'use strict';
define(function(require) {

   /**
    * Probability Distributions
    *
    * The correspondence between the possible values of a random
    * variable and their associated probabilities is referred to
    * as a *probability distribution*. Random variables can be *discrete*,
    * that is, taking any of a specified finite or countable list of values,
    * or *continuous*, taking any value in a specified subset
    * of the real numbers (typically, an interval).
    * This module provides some of the most commonly used discrete
    * and continuous probability distributions.
    *
    * If $X$ is a discrete random variable with sample space $S$
    * (the set of all possible values), then $X$ can be
    * specified by a probability mass function (pmf) $f:\,S \rightarrow [0,1]$.
    * The probability of $X$ taking the single value $x$ is given by $f(x)$:
    * $$\textrm{Pr}(X = x) = f(x)$$ and the total mass is 1:
    * $$\textrm{Pr}(X \in S) = \sum_{x \in S} f(x) = 1$$
    *
    * If $X$ is a continuous random variable, then $X$ can be specified
    * by a probability
    * density function (pdf) $f:\,\Re \rightarrow [0, \infty)$. The
    * probability of $X$
    * taking a value in a given interval $[a, b]$ is the corresponding
    * area beneath the density curve:
    * $$\textrm{Pr}(a \leq X \leq b) = \int_{a}^{b} f(x) \, dx$$
    *
    * The probability for a continuous random variable $X$ to take
    * any single value
    * is thus zero.
    * and the total mass is 1:
    * $$\textrm{Pr}(-\infty \leq X \leq \infty) = \int_{-\infty}^{\infty} f(x) \, dx = 1$$
    *
    * The cumulative distribution function (cdf) $F(x)$ is defined by
    * $F(x) = \textrm{Pr}(X \leq x)$; $F(x)$ is is often referred to as the
    * *left-tail probability* for the cutoff value $x$.
    * In the discrete case, where
    * $f(x)$ is the probability mass function and $S$ is
    * the sample space (set of all possible values for $X$),
    * the cdf is defined as a sum:
    * $$ F(x) = \sum_{t \in S \text{, } t \leq x}f(t)$$
    *
    * In the continuous case, where $f(x)$ is the probability density
    * function, the cdf is defined as an integral:
    * $$ F(x) = \int_{-\infty}^{x} f(t) \, dt$$
    *
    * For $0 < p < 1$, the *quantile function* $Q(p)$ returns the smallest
    * value of $X$ for which the left-tail probability is at least $p$:
    * $$Q(p) = \textrm{inf} \\{x : p \leq F(x) \\}$$
    *
    * When the cdf $F$ is continuous and strictly monotone,
    * $Q = F^{-1}$; that is,
    * $$Q(p) = x \text{ such that } \textrm{Pr}(X \leq x) = p$$
    * This explains the usage of *inverse cdf* as a synonym for *quantile
    * function*.
    *
    * If the sample space $S$ has a minimum, then
    * $Q(0) = \textrm{min}(S)$; similarly, if the sample space $S$ has a maximum, then
    * $Q(1) = \textrm{max}(S)$.
    *
    * A single observation from a random variable is referred to as a *variate*
    * or *deviate*.  Each of the distributions provided in this module includes
    * a method for generating random variates.
    *
    * @module distributions
    * @memberof panthrMath
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var distributions, mixin;

   mixin = require('./utils.js').mixin;

   distributions = mixin({},
      /** Constructor for finite distributions. See `module:finite`. */
      require('./distributions/finite'),
      /** The normal distribution. See `module:normal`. */
      require('./distributions/normal'),
      /** The gamma distribution. See `module:gamma`. */
      require('./distributions/gamma'),
      /** The chi squared distribution. See `module:chisq`. */
      require('./distributions/chisq'),
      /** The Poisson distribution. See `module:poisson`. */
      require('./distributions/poisson'),
      /** The t distribution. See `module:t`. */
      require('./distributions/t'),
      /** The uniform distribution. See `module:uniform`. */
      require('./distributions/uniform'),
      /** The beta distribution. See `module:beta`. */
      require('./distributions/beta'),
      /** The exponential distribution. See `module:exponential`. */
      require('./distributions/exp'),
      /** The binomial distribution. See `module:geometric`. */
      require('./distributions/geom'),
      /** The binomial distribution. See `module:binomial`. */
      require('./distributions/binom')
   );

   return distributions;

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
