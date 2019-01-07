package com.github.jonnylaw.gp

import breeze.linalg.DenseVector
import GaussianProcess._
import breeze.stats.distributions._
import breeze.linalg.{DenseVector, DenseMatrix, diag, cholesky}

sealed trait MeanParameters {
  def map(f: Double => Double): MeanParameters

  def toList: List[Double]

  def toMap: Map[String, Double]

  def add(that: MeanParameters): MeanParameters
}

case class Plane(beta: DenseVector[Double]) extends MeanParameters {
  def map(f: Double => Double): MeanParameters =
    MeanParameters.plane(beta mapValues f)

  def add(that: MeanParameters) = that match {
    case Plane(beta1) => MeanParameters.plane(beta + beta1)
    case _            => throw new Exception("Can't add different shaped parameters")
  }

  def toList = beta.data.toList

  def toMap = beta.data.zipWithIndex.map { case (b, i) =>
    (s"beta_$i" -> b)
  }.toMap
}

case object Zero extends MeanParameters {
  def map(f: Double => Double): MeanParameters =
    MeanParameters.zero

  def add(that: MeanParameters) = that match {
    case Zero => MeanParameters.zero
    case _    => throw new Exception("Can't add different shaped parameters")
  }

  def toList = Nil

  def toMap = Map.empty[String, Double]
}

object MeanParameters {
  def plane(beta: DenseVector[Double]): MeanParameters = {
    Plane(beta)
  }

  def zero: MeanParameters = Zero

  /**
    * Make a design matrix from the location vectors
    */
  def makeDesignMatrix(x: Vector[Location[Double]]): DenseMatrix[Double] = {
    val xs: Vector[Vector[Double]] = x.map(_.toVector)

    val n = x.size
    val m = xs.head.size
    val res = new DenseMatrix(n, m, xs.flatten.toArray)

    DenseMatrix.horzcat(DenseVector.ones[Double](res.rows).toDenseMatrix.t, res)
  }

  /**
    * Sample the value of the mean parameters from the full conditional distribution
    * using a conditionally conjugate prior in a Gibbs Step
    * @param prior Gaussian the conditionally conjugate prior
    * @param obs a collection of multivariate observations
    */
  def samplePlane(
    prior: Gaussian,
    obs:   Vector[Data],
    xs:    Vector[Location[Double]],
    dist:  (Location[Double], Location[Double]) => Double,
    p:     Parameters): MeanParameters = {

    val covFn = KernelFunction.apply(p.kernelParameters)
    // build the design matrix
    val x = makeDesignMatrix(xs)
    val nugget = diag(DenseVector.fill(xs.size)(1e-6))
    val kxx = KernelFunction.buildCov(xs, covFn, dist) + nugget
    val l = cholesky(kxx)

    // calculate the precision and cholesky factor
    val priorPrec = diag(DenseVector.fill(x.cols)(1.0 / prior.variance))
    val prec = x.t * (kxx \ x) + priorPrec

    val y = DenseVector(obs.map(_.y).toArray) 

    val u: DenseVector[Double] = Predict.forwardSolve(l, y)
    val mean: DenseVector[Double] = prec \ (priorPrec * DenseVector.fill(
      x.cols)(prior.mean) + (x.t * u))

    val root = cholesky(prec)
    val z = DenseVector.fill(mean.size)(Gaussian(0, 1).draw)
    val beta = mean + root \ z

    MeanParameters.plane(beta)
  }

  /**
    * Sample mean hyperparameters
    * @param delta the variance of the Gaussian proposal distribution
    * @param obs a collection of multivariate observations
    * @param prior the prior distribution of the hyperparameters
    */
  def sample(
    obs:   Vector[Data],
    xs:    Vector[Location[Double]],
    dist:  (Location[Double], Location[Double]) => Double,
    prior: Gaussian
  ) = { p: Parameters => p.meanParameters match {
    case Plane(beta) =>
      Rand.always(p.copy(meanParameters = samplePlane(prior, obs, xs, dist, p)))
    case Zero =>
      Rand.always(p)
  }}

  def vectorToParams(p: MeanParameters, pActual: Vector[Double]) =
    p match {
      case Plane(beta) =>
        plane(DenseVector(pActual.toArray))
      case Zero => MeanParameters.zero
    }
}
