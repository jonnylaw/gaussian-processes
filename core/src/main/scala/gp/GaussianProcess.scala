package com.github.jonnylaw.gp

import breeze.stats.distributions._
import breeze.linalg.{DenseVector, DenseMatrix, cholesky, sum, diag, eigSym}
import breeze.numerics.log
import cats.implicits._

object GaussianProcess {

  /**
    * Parameters of a Gaussian Process model
    * @param sigma the noise standard deviation of observations of the process
    * @param kernelParameters the parameters of the covariance function
    */
  case class Parameters(
      meanParameters: MeanParameters,
      kernelParameters: Vector[KernelParameters]
  ) {

    def map(f: Double => Double): Parameters = {
      Parameters(meanParameters map f, kernelParameters.map(_.map(f)))
    }

    def toList: List[Double] =
      meanParameters.toList ++ kernelParameters.toList.flatMap(_.toList)
  }

  /**
    * A single measurement of a process in multiple dimensions
    * @param x a location (eg. a Latitude and Longitude pair)
    * @param y a measurement
    */
  case class Data(x: Location[Double], y: Double)

  /**
    * A realisation of a Gaussian Process at
    */
  case class TimePoint(t: Double, ys: Vector[Data])

  /**
    * Calculate the pairwise distance matrix for a vector of positions
    * @param xs a vector of locations
    * @param dist a distance metric
    * @return a symmetrix matrix containing the pairwise distances between all
    * locations in xs
    */
  def distanceMatrix(xs: Vector[Location[Double]],
                     dist: (Location[Double], Location[Double]) => Double) = {

    val n = xs.length

    val res: Vector[Double] = for {
      x <- xs
      y <- xs
    } yield dist(x, y)

    new DenseMatrix(n, n, res.toArray)
  }

  /**
    * Given a covariance and mean function make
    * a draw from the Gaussian Process prior at specified locations
    */
  def draw(xs: Vector[Location[Double]],
           dist: (Location[Double], Location[Double]) => Double,
           p: Parameters)(implicit rand: RandBasis = Rand) = {

    val nugget = diag(DenseVector.fill(xs.size)(1e-3))
    val covariance = KernelFunction.buildCov(
      xs,
      KernelFunction.apply(p.kernelParameters),
      dist) + nugget

    val mu = DenseVector(xs.map(MeanFunction.apply(p.meanParameters)).toArray)

    val root = eigSym(covariance)

    val x = DenseVector.rand(mu.length, rand.gaussian(0, 1))
    mu + (root.eigenvectors * diag(root.eigenvalues.mapValues(math.sqrt)) * x)
  }

  def vecToData(vec: DenseVector[Double], xs: Vector[Location[Double]]) = {
    (vec.data.toVector zip xs) map { case (f, x) => Data(x, f) }
  }

  def drawData(
      xs: Vector[Location[Double]],
      dist: (Location[Double], Location[Double]) => Double,
      p: Parameters
  ) = {

    val covariance = KernelFunction.buildCov(
      xs,
      KernelFunction.apply(p.kernelParameters),
      dist)
    val mean = DenseVector(xs.map(MeanFunction.apply(p.meanParameters)).toArray)

    for {
      fx <- MultivariateGaussian(mean, covariance)
    } yield vecToData(fx, xs)
  }

  /**
    * Calculate the marginal log-likelihood of a single observation
    * of a Gaussian process
    */
  def loglikelihood(observed: Vector[Data],
                    dist: (Location[Double], Location[Double]) => Double) = {
    p: Parameters =>
    val covFn = KernelFunction(p.kernelParameters)
    val xs = observed.map(_.x)
    val n = xs.size

    // covariance of observed
    val nugget = diag(DenseVector.fill(xs.size)(1e-3))
    val kxx = KernelFunction.buildCov(xs, covFn, dist) + nugget

    val meanFn = MeanFunction.apply(p.meanParameters)
    val ys = DenseVector(observed.map(_.y).toArray) - DenseVector(
      xs.map(meanFn).toArray)
    val l = cholesky(kxx)
    val u = Predict.forwardSolve(l, ys)

    -0.5 * u dot u - sum(log(diag(l))) - n * 0.5 * log(2 * math.Pi)
  }

  /**
    * Take a random sampling of locations between two points in one-dimension
    */
  def samplePoints(start: Double, end: Double, n: Int): Vector[Double] = {
    Uniform(start, end).sample(n).sorted.toVector
  }

  /**
    * Draw from a conditional Gaussian distribution in a more efficient way
    * because we only have to calculate the cholesky of the covariance of the
    * prior distribution which only has to be computed once for every new sample X | Y
    * http://www.stats.ox.ac.uk/~doucet/doucet_simulationconditionalgaussian.pdf
    * Sampling X | Y
    * This can be used when we require several draws from the same posterior
    * @param xs a vector of locations at which to determine f(x)
    * @param ys a vector of noisy observations of y = f(x) + eps
    */
  def efficientDraw(
      xs: Vector[Location[Double]],
      ys: Vector[Data],
      dist: (Location[Double], Location[Double]) => Double,
      p: Parameters,
  )(priorSample: Rand[DenseVector[Double]]) = {

    val covFn = KernelFunction.apply(p.kernelParameters)
    val ysx = ys.map(_.x)
    val kyy = KernelFunction.buildCov(ysx, covFn, dist)
    val kxy = KernelFunction.buildDistCov(xs, ysx, covFn, dist)

    for {
      z <- priorSample
      x = z.slice(0, xs.size)
      y = z.slice(xs.size, -1)
      yvec = DenseVector(ys.map(_.y).toArray)
    } yield x + (kyy \ kxy) * (y - yvec)
  }

  def mllGradient(
    observed: Vector[Data],
    dist: (Location[Double], Location[Double]) => Double
  )(p: Parameters): Array[Double] = {

    val covFn = KernelFunction(p.kernelParameters)
    val xs = observed.map(_.x)

    // covariance of observed
    val nugget = diag(DenseVector.fill(xs.size)(1e-3))
    val kxx = KernelFunction.buildCov(xs, covFn, dist) + nugget

    val meanFn = MeanFunction.apply(p.meanParameters)
    val ys = DenseVector(observed.map(_.y).toArray) - DenseVector(
      xs.map(meanFn).toArray)

    val grad = KernelParameters.tangentMatrix(observed, dist, p.kernelParameters)

    val alpha = kxx \ ys

    grad.map { g => 0.5 * sum(diag(alpha * alpha.t * g - kxx \ g)) }.toArray
  }
}
