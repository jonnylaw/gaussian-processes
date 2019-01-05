---
layout: docs
title: "Introduction to the Gaussian Process"
---

A Gaussian Process can be used to approximate a non-linear function $y = f(x)$. To get started with the Gaussian Process package, first import the library:

```tut:silent
import com.github.jonnylaw.gp._
```

Finitely many points of a Gaussian process are distributed according to a
multivariate Normal (MVN) distribution. The MVN is parameterised by a mean and
[covariance function](covariance_function.html). A covariance function is
determined by the distance between points and some static parameters. The
distance function can be any suitable function from a pair of locations to a
distance measurement `(Location[Double], Location[Double]) => Double`. Euclidean
distance is a suitable distance function when dealing with small scale spatial data:

$$\textrm{d}(x, y) = \sqrt{x^2 - y^2}$$

Other suitable distance functions include Manhattan distance or Great Circle
distance for measurements on global scale.

In order to simulate data from a Gaussian Process, first specify a vector of locations:

```tut
val xs = GaussianProcess.samplePoints(-10.0, 10.0, 300).
  map(One.apply)
```

The function `samplePoints` uniformly samples a vector of 300 Doubles between -10
and 10. The sampled points are then put into a `Location` object `One`. Next a
matrix representing the pairwise distances between each point can be calculated:

```tut
val m = GaussianProcess.distanceMatrix(xs, Location.euclidean)
```

Next a covariance function can be selected, for instance the squared exponential
covariance function:

$$K(x, y) = h\exp\left\{ \frac{1}{\sigma^2}\textrm{d}(x, y)^2 \right\}$$

$\sigma$ is known as the length scale, large values of sigma indicate that large
changes in distance are required for the covariance to change significantly. $h$
controls the amount of change between two points. The squared exponential
function can be applied to the distance matrix using a `map`

```tut
val covFn: Double => Double = KernelFunction.squaredExponential(h = 3.0, sigma =
5.0)

val covMat = m.map(covFn)
```

Then a draw from a zero-mean Gaussian Process prior with covariance matrix
`covMat` at locations `xs` can be performed:

```tut
val zero = DenseVector.zero[Double](300)
val ys = breeze.stats.distributions.MultivariateGaussian(zero, covMat).draw
```

Then the data can be plotted:

```tut
Plot.scatterPlot(GaussianProcess.vecToData(xs, ys))
```
