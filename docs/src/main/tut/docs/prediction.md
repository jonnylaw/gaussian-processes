---
layout: docs
title: "Fitting a Gaussian Process"
---

Given values of the hyper parameters of the covariance and mean functions. Then
we can determine the posterior distribution of the function $f(x)$ given paired
data $(x, y)$ where $x$ is in the domain of the function and $y$ is in the range
of the function. Typically, $x$ will be a location or time and $y$ is a process
we want to learn about.

For a continuous outcome, $y$ the likelihood and prior is Gaussian, resulting in
a posterior distribution which is also Gaussian and can be derived exactly.
Given a location with an unknown measurement, $y(x^\star)$, then the posterior
is Gaussian and can be determined analytically by calculating the mean and variance:

\begin{align*}
\mathbb{E}(y(x^\star)) &= K(x_\star, \textbf{x})^T(K(\textbf{x}, \textbf{x}) + I_n\sigma_y)^{-1}\textbf{y} \\
\textrm{Var}(y(x^\star)) &= K(x_\star, x_\star) - K(x_\star, \textbf{x})^T(K(\textbf{x}, \textbf{x}) + I_n\sigma_y^2)^{-1}K(x_\star, \textbf{x})
\end{align*}

Suppose that we want to determine the posterior distribution of the function
$f(x)$ by observing finitely many points from the simulation in the [GP
introduction](introduction.html). This can be simulated by specifying the
parameters and a distance function:

```tut
import com.github.jonnylaw.gp._

val params = GaussianProcess.Parameters(
  MeanParameters.zero,
  Vector(KernelParameters.se(3.0, 5.5))
)
val dist = Location.euclidean _

val xs = GaussianProcess.samplePoints(-10.0, 10.0, 300).map(One.apply)
val ys = GaussianProcess.draw(xs, dist, params)
```

Then observing only every 15th point:

```tut
val observed = GaussianProcess.vecToData(ys, xs).
    zipWithIndex.
    filter { case (_, i) => (i + 1) % 15 == 0 }.
    map(_._1)
```

In order to determine the posterior distribution of the function, we specify a
list of test-points which are in the same domain as the function and calculate
the posterior mean and variance given the observed data:

```tut
implicit val integralD = scala.math.Numeric.DoubleAsIfIntegral
val testPoints = Vector.range(-10.0, 10.0, 0.01).map(One(_))
// draw from the GP posterior at values of testPoints
val fitted = Predict.fit(testPoints ++ observed.map(_.x), observed, dist, params)
```

This calculates the mean and covariance at each of the test points. This can be
used to determine posterior probability intervals for the function which can then be plotted:

```tut
import com.cibo.evilplot.plot.aesthetics.DefaultTheme._

Plot.gpPlot(fitted)
com.cibo.evilplot.plot.Overlay(Plot.scatterPlot(observed)).
  render()
//  write(new java.io.File("figures/fitted_gp.png"))
```

![Fitted GP](figures/fitted_gp.png)

