package gp.core

import breeze.stats.distributions._
import breeze.linalg.{DenseVector, DenseMatrix, cholesky, diag}
import GaussianProcess._
// import com.github.fommil.netlib.LAPACK.{getInstance => lapack}
import com.github.fommil.netlib.BLAS.{getInstance => blas}

object Predict {
  /**
    * Helper function for fit chol
    * @param newx a new location
    * @param xs locations with measurements
    * @param covFn a covariance function from Double => Double
    * @param dist a distance function from a pair of locations to double
    * @return a vector containing the value of the covariance function
    * given the distance from each existing location to the new location
    */
  def buildDistVec(
    newx:  Location[Double],
    xs:    Vector[Location[Double]],
    covFn: Double => Double,
    dist:  (Location[Double], Location[Double]) => Double
  ): DenseVector[Double] = {

    DenseVector(xs.map(d => covFn(dist(d, newx))).toArray)
  }

  /**
    * Forward solve a lower-triangular linear system
    * with a single RHS
    *
    * @param A A lower-triangular matrix
    * @param y A single vector RHS
    *
    * @return The solution, x, of the linear system A x = y
    */
  def forwardSolve(
    A: DenseMatrix[Double],
    y: DenseVector[Double]): DenseVector[Double] = {
    val yc = y.copy
    blas.dtrsv("L", "N", "N", A.cols, A.toArray, A.rows, yc.data, 1)
    yc
  }

  /**
    * Algorithm 2.1 from GPML
    * An efficient way of learning a Gaussian Process at finitely many
    * points newXs
    * @param newXs the new locations we are interesting in determining the
    * value of the posterior of f
    * @param observed the currently observed (noisy) values of the function
    * @param mod a Gaussian Process Model
    * @param p Gaussian process parameters
    */
  def fit(
    newXs:    Vector[Location[Double]],
    observed: Vector[Data],
    dist:     (Location[Double], Location[Double]) => Double,
    p:        Parameters): Vector[(Location[Double], Gaussian)] = {

    val covFn = KernelFunction.apply(p.kernelParameters)

    val x = observed.map(_.x)
    // add a small amount for numerical stability (tikhonov regularization)
    val nugget = diag(DenseVector.fill(x.size)(1e-6))
    val kxx = KernelFunction.buildCov(x, covFn, dist) + nugget
    val l = cholesky(kxx)

    for {
      xs <- newXs

      // calculate covariances, between the new points and existing
      kxy = buildDistVec(xs, x, covFn, dist)

      // between the point and itself
      kyy = covFn(dist(xs, xs)) //

      ys = DenseVector(observed.map(_.y).toArray)

      // calculate the mean and covariance
      // alpha = l.t \ (l \ ys)
      // v = l \ kxy
      // mean = kxy.t * alpha
      // cov = kyy - v.t * v

      u = forwardSolve(l, kxy)
      v = forwardSolve(l, ys)
      mean = u.t * v
      cov = kyy - u.t * u
    } yield (xs, Gaussian(mean, math.sqrt(cov)))
  }

  def predict(
    fitted:   Vector[Gaussian],
    interval: Double) = {

    fitted.map(f =>
      (f.mean, Summarise.getInterval(f.mean, f.variance, 1 - interval),
        Summarise.getInterval(f.mean, f.variance, interval)))
  }

  def buildData(
    xs: Vector[Location[Double]],
    ys: Vector[Double]): Vector[Data] = {

    (xs zip ys) map { case (x, y) => Data(x, y) }
  }
}
