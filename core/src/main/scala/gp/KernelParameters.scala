package gp.core

import breeze.stats.distributions._
import breeze.linalg.{DenseMatrix, DenseVector}

import GaussianProcess._

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
  def paramsToDenseVector(p: Parameters) =
    DenseVector(p.toList.toArray)

  /**
    * Transform a dense vector to parameters
    */
  def vectorToParams(p: Parameters, pActual: DenseVector[Double]) = {
    val n = p.meanParameters.toList.size
    val ps = pActual.data.toVector
    Parameters(
      MeanParameters.vectorToParams(p.meanParameters, ps.take(n)),
      KernelParameters.vectorToParams(p.kernelParameters, ps.drop(n))
    )
  }

  /**
    * Sample GP parameters using HMC
    * @param ys a vector of observations of a GP
    */
  def sampleHmc(
    ys:   Vector[Data],
    dist: (Location[Double], Location[Double]) => Double,
    init: Parameters,
    ll:   Parameters => Double,
    m:    Int,
    l:    Int,
    eps:  Double) = {

    def newLl(p: DenseVector[Double]) = {
      val params = vectorToParams(init, p)
      val const = PernelParameters.
      ll(params.copy(kernelParameters = kp))
    }
    def newGrad(p: DenseVector[Double]) = {
      val params = vectorToParams(init, p)

      mllGradient(ys, dist)(params.copy(kernelParameters = kp)) + 
    }

    val kp = KernelParameters.unconstrainParams(init.kernelParameters)
    val theta = paramsToDenseVector(init.copy(kernelParameters = kp))
    val initState = HmcState(0, theta, 0)
    MarkovChain(initState)(Hmc(m, l, eps, newGrad, newLl).step)
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
      val h = Hmc.boundedBelow(0)(hu)
      val sc = Hmc.boundedBelow(0)(su)
      val s = sc.constrained

      Vector(math.exp(-(dist * dist) / (s * s)) + h.derivative,
        - (2 * h.constrained / (s * s * s)) * math.exp(-(dist * dist) / (s * s)) + sc.derivative)

    case White(su) => {
      val s = Hmc.boundedBelow(0)(su)
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
}
