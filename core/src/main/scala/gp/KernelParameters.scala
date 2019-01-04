package com.github.jonnylaw.gp

import breeze.stats.distributions._
import breeze.linalg.{DenseMatrix, DenseVector}

import GaussianProcess._
import breeze.linalg.DenseMatrix
import cats.implicits._
import math._

sealed trait KernelParameters { self =>
  def map(f: Double => Double): KernelParameters

  def toList: List[Double]

  def add(that: KernelParameters): KernelParameters
}

case class SquaredExp(h: Double, sigma: Double) extends KernelParameters {
  def map(f: Double => Double): KernelParameters = {
    KernelParameters.se(f(h), f(sigma))
  }

  def toList = h :: sigma :: Nil

  def add(that: KernelParameters) = that match {
    case SquaredExp(h1, sigma1) =>
      KernelParameters.se(h + h1, sigma + sigma1)
    case _ => throw new Exception
  }
}

case class Matern(sigma: Double, nu: Double, l: Double) extends KernelParameters {
  def map(f: Double => Double): KernelParameters = {
    KernelParameters.matern(f(sigma), f(nu), f(l))
  }

  def toList = sigma :: nu :: l :: Nil

  def add(that: KernelParameters) = that match {
    case Matern(sigma1, nu1, l1) =>
      KernelParameters.matern(sigma + sigma1, nu + nu1, l + l1)
    case _ => throw new Exception
  }
}

case class White(sigma: Double) extends KernelParameters {
  def map(f: Double => Double): KernelParameters = {
    KernelParameters.white(f(sigma))
  }

  def toList = sigma :: Nil

  def add(that: KernelParameters) = that match {
    case White(sigma1) =>
      KernelParameters.white(sigma + sigma1)
    case _ => throw new Exception
  }
}

object KernelParameters {
  /**
    * Smart constructor for the squared exponential covariance parameters
    */
  def se(h: Double, sigma: Double): KernelParameters = {
    SquaredExp(h, sigma)
  }

  def matern(sigma: Double, nu: Double, l: Double): KernelParameters = {
    Matern(sigma, nu, l)
  }

  def white(sigma: Double): KernelParameters = {
    White(sigma)
  }

  /**
    * Transform parameters to a dense vector
    */
  def paramsToArray(p: Parameters) =
    p.toList.toArray

  /**
    * Transform a dense vector to parameters
    */
  def arrayToParams(p: Parameters, pActual: Array[Double]): Parameters = {
    val n = p.meanParameters.toList.size
    val ps = pActual.toVector
    Parameters(
      MeanParameters.vectorToParams(p.meanParameters, ps.take(n)),
      KernelParameters.vectorToParams(p.kernelParameters, ps.drop(n))
    )
  }

  case class HmcState(
    iter:     Int,
    theta:    DenseVector[Double],
    accepted: Int)

  /**
    * Sample GP parameters using HMC
    * @param ys a vector of observations of a GP
    */
  def sampleHmc(
    ys:   Vector[Data],
    dist: (Location[Double], Location[Double]) => Double,
    init: Parameters,
    ll:   Parameters => Double,
    m:    DenseMatrix[Double],
    l:    Int,
    eps:  Double): Process[Array[Double]] = {

    def newLl(p: Array[Double]) = {
      val params = arrayToParams(init, p)
      val kp = KernelParameters.constrainParams(params.kernelParameters)
      ll(params.copy(kernelParameters = kp))
    }
    def newGrad(p: Array[Double]) = {
      val params = arrayToParams(init, p)
      mllGradient(ys, dist)(params)
    }

    // unconstrain the parameters to they are on the whole real line
    val kp = KernelParameters.unconstrainParams(init.kernelParameters)
    val theta = paramsToArray(init.copy(kernelParameters = kp))
    MarkovChain(theta) { theta =>
      for {
        (tp, _) <- Hmc.step(eps, l, m, newLl, newGrad)(theta)
      } yield tp
    }
  }

  /**
    * Perform sampling using empirical HMC, requiring no tuning of step
    * size or leapfrog steps
    */
  def sampleEhmc(
    ys: Vector[Data],
    dist: (Location[Double], Location[Double]) => Double,
    theta: Parameters,
    ll: Parameters => Double,
    l0: Int,
    mu: Double,
    warmupIterations: Int) = {

    def newLl(p: Array[Double]) = {
      val params = arrayToParams(theta, p)
      val kp = KernelParameters.constrainParams(params.kernelParameters)
      ll(params.copy(kernelParameters = kp))
    }
    def newGrad(p: Array[Double]) = {
      val params = arrayToParams(theta, p)
      mllGradient(ys, dist)(params)
    }

    // unconstrain the parameters to they are on the whole real line
    val kp = KernelParameters.unconstrainParams(theta.kernelParameters)
    val initTheta = paramsToArray(theta.copy(kernelParameters = kp))
    val m = DenseMatrix.eye[Double](initTheta.size)
    Ehmc.sample(l0, mu, m, newLl, newGrad, warmupIterations, initTheta)
  }


  /**
    * Sample the observation variance from a Gamma distribution
    * Gibbs step using a conditionally-conjugate distribution
    * @param prior the conditional conjugate prior distribtion for the
    * measurement noise variance
    * @param ys a vector of data points
    * @param fx the currently sampled value of the function state
    * @return a distribution over the measurement noise
    */
  def samplePrecY(
    prior: Gamma,
    ys:    Vector[Data],
    fx:    Vector[Data]) = {

    val ssy = (ys zip fx).map {
      case (y, f) => (y.y - f.y) * (y.y - f.y)
    }.sum

    val shape = prior.shape + ys.size * 0.5
    val scale = prior.scale + 0.5 * ssy

    Gamma(shape, scale)
  }

  /**
    * Sample kernel hyperparameters using the metropolis-hastings algorithm
    * @param obs the observations
    * @param dist a distance function
    * @param prior the prior distribution of the hyperparameters
    * @param prop a symmetric proposal distribution
    */
  def sample(
    obs:   Vector[Data],
    dist:  (Location[Double], Location[Double]) => Double,
    prior: Vector[KernelParameters] => Double,
    prop:  Vector[KernelParameters] => Rand[Vector[KernelParameters]]) = {

    val proposal = (p: Parameters) =>
    for {
      kp <- prop(p.kernelParameters)
    } yield p.copy(kernelParameters = kp)

    val ll = (p: Parameters) => prior(p.kernelParameters) +
      GaussianProcess.loglikelihood(obs, dist)(p)

    MarkovChain.Kernels.metropolis(proposal)(ll)
  }

  /**
    * Transform the parameters to the whole of the real line
    */
  def unconstrainParams(p: Vector[KernelParameters]): Vector[KernelParameters] = p map {
    case SquaredExp(h, s) => se(log(h), log(s))
    case White(s) => white(log(s))
    case Matern(_, _, _) => throw new Exception("Not implemented yet")
  }

  /**
    * Restrict the Kernel parameters to an appropriate domain
    */
  def constrainParams(p: Vector[KernelParameters]): Vector[KernelParameters] = p map {
    case SquaredExp(h, s) => se(exp(h), exp(s))
    case White(s) => white(exp(s))
    case Matern(_, _, _) => throw new Exception("Not implemented yet")
  }

  /**
    * Calculate the gradient of the unconstrained Kernel Parameters
    * for a given distance
    */
  def gradient(p: Vector[KernelParameters])(dist: Double): Vector[Double] = p flatMap {
    case SquaredExp(hu, su) =>
      val h = boundedBelow(0)(hu)
      val sc = boundedBelow(0)(su)
      val s = sc.constrained

      Vector(math.exp(-(dist * dist) / (s * s)) + h.derivative,
        - (2 * h.constrained / (s * s * s)) * math.exp(-(dist * dist) / (s * s)) + sc.derivative)

    case White(su) => {
      val s = boundedBelow(0)(su)
      Vector(1.0 + s.derivative)
    }
    case Matern(_, _, _) => throw new Exception("Not implemented yet")
  }

  /**
    * Build the tangent matrices (a matrix of partial derivatives)
    */
  def tangentMatrix(
    observed: Vector[Data],
    dist: (Location[Double], Location[Double]) => Double,
    p:    Vector[KernelParameters]): Vector[DenseMatrix[Double]] = {

    val xs = observed.map(_.x)

    val m = distanceMatrix(xs, dist)
    val res = m.data.toVector map (gradient(p))

    res.sequence.map (ms => new DenseMatrix(m.rows, m.cols, ms.toArray))
  }

  /**
    * Transform a DenseVector containing the values of the parameter
    * into a set of Kernel parameters
    */
  def vectorToParams(p: Vector[KernelParameters],
                     pActual: Vector[Double]): Vector[KernelParameters] = {

    p.foldLeft((Vector[KernelParameters](), pActual)) {
        case ((acc, pa), ps) =>
          ps match {
            case _: SquaredExp =>
              (acc :+ se(pa.head, pa(1)), pa.drop(2))
            case _: White =>
              (acc :+ white(pa.head), pa.tail)
            case _: Matern =>
              (acc :+ matern(pa.head, pa(1), pa(2)), pa.drop(3))
          }
      }
      ._1
  }

  def logistic(x: Double): Double =
    1.0 / (1.0 + exp(-x))

  def logit(p: Double): Double =
    log(p / (1 - p))

  def softplus(x: Double): Double =
    log1p(exp(x))

  /**
    * Parameters in an HMC algorithm
    */
  case class Parameter(
    unconstrained: Double,
    constrained: Double,
    logJacobian: Double,
    derivative: Double,
    unconstrain: Double => Double,
    constrain: Double => Double
  )

  def unbounded(x: Double) =
    Parameter(x, x, 0, 0, x => x, x => x)

  def bounded(min: Double, max: Double)(x: Double) = Parameter(
    x,
    logistic(x) * (max - min) + min,
    log(max - min) - x + 2.0 * log(logistic(x)),
    -1.0 + 2 * exp(-x) / (1.0 + exp(-x)),
    y => logit((y - min) / (max - min)),
    (y) => logistic(y) * (max - min) + min,
    )

  def boundedBelow(min: Double)(x: Double) = Parameter(
    x,
    exp(x) + min,
    x,
    1.0,
    y => log(y - min),
    y => exp(y) + min)

  def boundedAbove(max: Double)(x: Double) = Parameter(
    x,
    max - exp(-x),
    -x,
    -1.0,
    y => -log(max - y),
    y => max - exp(-y))

}

