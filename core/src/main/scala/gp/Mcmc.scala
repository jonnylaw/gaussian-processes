package gp.core

import GaussianProcess._
import breeze.stats.distributions._
import breeze.linalg.{DenseVector, DenseMatrix, cholesky, diag, inv}
import cats.implicits._
import com.github.fommil.netlib.BLAS.{getInstance => blas}
import cats.data.Kleisli

/**
  * Use Gibbs Sampling and Metropolis-Hastings
  * to determine the state and parameters of a model
  * group by locations in data
  */
object Mcmc {

  /**
    * The state of the MCMC sampler, containing the currently sampled
    * values of the parameters and values of the function
    * @param p the current parameters
    * @param fx realisations of the function we are learning
    */
  case class State(p: Parameters, fx: Vector[Data])

  /**
    * Backsolve an upper-triangular linear system
    * with a single RHS
    *
    * @param A An upper-triangular matrix
    * @param y A single vector RHS
    *
    * @return The solution, x, of the linear system A x = y
    */
  def backSolve(A: DenseMatrix[Double],
                y: DenseVector[Double]): DenseVector[Double] = {
    val yc = y.copy
    blas.dtrsv("U", "N", "N", A.cols, A.toArray, A.rows, yc.data, 1)
    yc
  }

  /**
    * Backsolve an upper-triangular linear system
    * with multiple RHSs
    *
    * @param A An upper-triangular matrix
    * @param Y A matrix with columns corresponding to RHSs
    *
    * @return Matrix of solutions, X, to the linear system A X = Y
    */
  def backSolve(A: DenseMatrix[Double],
                Y: DenseMatrix[Double]): DenseMatrix[Double] = {
    val yc = Y.copy
    blas.dtrsm("L",
               "U",
               "N",
               "N",
               yc.rows,
               yc.cols,
               1.0,
               A.toArray,
               A.rows,
               yc.data,
               yc.rows)
    yc
  }

  /**
    * Sample the Gaussian Process from the prior, this will result in slow mixing
    */
  def sampleStatePrior(xs: Vector[Location[Double]],
                       dist: (Location[Double], Location[Double]) => Double) = {
    s: State =>
      GaussianProcess.drawData(xs, dist, s.p) map (newfx => s.copy(fx = newfx))
  }

  /**
    *
    */
  def step(ys: Vector[Data],
           xs: Vector[Location[Double]],
           priorKern: Vector[KernelParameters] => Double,
           priorMean: Gaussian,
           propKern: Vector[KernelParameters] => Rand[Vector[KernelParameters]],
           dist: (Location[Double], Location[Double]) => Double) = {

    Kleisli(KernelParameters.sample(ys, dist, priorKern, propKern)) andThen
      Kleisli(MeanParameters.sample(ys, xs, dist, priorMean))
  }

  /**
    * Sample the Parameters
    */
  def sample(
      ys: Vector[Data],
      priorKern: Vector[KernelParameters] => Double,
      priorMean: Gaussian,
      propKern: Vector[KernelParameters] => Rand[Vector[KernelParameters]],
      dist: (Location[Double], Location[Double]) => Double,
      initP: Parameters
  ) = {

    val xs = ys.map(_.x).distinct

    MarkovChain(initP)(step(ys, xs, priorKern, priorMean, propKern, dist).run)
  }
}
