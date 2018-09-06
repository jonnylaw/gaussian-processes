package gp.core

import breeze.numerics._
import breeze.linalg.DenseMatrix
import GaussianProcess._
import org.apache.commons.math3.special.Gamma

object KernelFunction {
  /**
    * Evaluate the squared exponential covariance kernel
    */
  def squaredExponential(h: Double, sigma: Double)(dist: Double): Double = {
    h * exp(- (dist * dist) / (sigma * sigma) )
  }

  /**
    * The Matern covariance function 
    * @param nu 
    * @param l
    */
  def matern(
    sigma: Double,
    nu:    Double,
    l:     Double)
    (dist: Double): Double = {

    val vlr = (sqrt(2 * nu) * dist) / l

    sigma * 1.0 / (math.pow(2, nu - 1) * Gamma.gamma(nu)) * math.pow(vlr, nu) * Bessel.i0(vlr)
  }

  /**
    * A white noise kernel representing measurement noise
    */
  def white(sigma: Double)(dist: Double): Double = {
    sigma
  }

  /**
    * Apply the Kernel function as a sum 
    * 
    * TODO: This could have other ways of combining the covariance functions
    * (convolution, product etc.)
    * 
    * @param ps a list of covariance function parameters
    * @return a function from distance to the evaluation of a sum covariance function
    */
  def apply(ps: Vector[KernelParameters]): Double => Double = { x => 
    ps.map(p => p match {
      case SquaredExp(h, s) => squaredExponential(h, s)(x)
      case Matern(s, nu, l) => matern(s, nu, l)(x)
      case White(s) => white(s)(x)
    }).sum
  }

  /**
    * Build a covariance matrix S = ((k_xx k_xy), (k_xy.t k_yy)
    */
  def buildCovMatrix(
    kxx: DenseMatrix[Double],
    kyy: DenseMatrix[Double],
    kxy: DenseMatrix[Double]): DenseMatrix[Double] = {

    val n = kxx.rows + kyy.rows
    val m = kxx.rows

    DenseMatrix.tabulate(n, n){ case (i, j) =>
      if (i < m & j < m) {
        kxx(i,j)
      } else if (i > m & j < m) {
        kxy(i-m, j)
      } else if (i < m & j > m) {
        kxy(i, j-m)
      } else {
        kyy(i-m, j-m)
      }
    }
  }

  /**
    * Build a covariance matrix of paired distances
    */
  def buildCov(
    xs:    Vector[Location[Double]],
    covFn: Double => Double,
    dist:  (Location[Double], Location[Double]) => Double
  ) = {
    distanceMatrix(xs, dist).map(covFn)
  }

  /**
    * Build distance from new points to existing observations
    */
  def buildDistCov(
    newXs: Vector[Location[Double]],
    xs:    Vector[Location[Double]],
    covFn: Double => Double,
    dist:  (Location[Double], Location[Double]) => Double) = {

    val n = newXs.size
    val m = xs.size
    val mat = DenseMatrix.zeros[Double](n, m)
    
    for {
      i <- 1 until n
      j <- 1 until m
    } mat(i, j) = covFn(dist(newXs(i), xs(j)))

    mat
  }
}
