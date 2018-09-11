package gp.core

import breeze.stats.distributions._
import breeze.linalg.{DenseVector, DenseMatrix}
import GaussianProcess._
import math._

/**
  * Discrete uniform distribution
  */
case class DiscreteUniform(min: Int, max: Int) extends Rand[Int] {
  def draw =
    min + scala.util.Random.nextInt(max - min + 1)
}

/**
  * TODO: Change leapfrog steps to handle constraints using section 5.1 of 
  * handbook of MCMC
  */
case class HmcState(
  theta:    DenseVector[Double],
  accepted: Int)

/**
  * Hamiltonian Monte Carlo 
  * prior distribution for the momentum parameters
  */
case class Hmc(
  m: Int,
  l: Int,
  eps: Double,
  gradient: DenseVector[Double] => DenseVector[Double],
  ll: DenseVector[Double] => Double) {

  /**
    * Prior distribution for the momentum variables
    */
  def priorPhi = {
    val zero = DenseVector.zeros[Double](m)
    MultivariateGaussian(zero, DenseMatrix.eye[Double](m))
  }

  /**
    * Update the value of the momentum
    */
  private def leapfrogHalfStep(
    theta: DenseVector[Double],
    phi: DenseVector[Double],
    eps: Double) = 
      phi + eps * 0.5 * gradient(theta)

  /**
    * Perform a leapfrog update with 1 step
    */
  def leapfrog(
    theta: DenseVector[Double],
    phi: DenseVector[Double],
    eps: Double) = {

    val p1 = leapfrogHalfStep(theta, phi, eps)
    val t1 = theta + eps * p1
    val p2 = leapfrogHalfStep(t1, p1, eps)

    (t1, p2)
  }

  private def leapfrogs(
    theta: DenseVector[Double],
    phi: DenseVector[Double],
    eps: Double): Stream[(DenseVector[Double], DenseVector[Double])] = {

    Stream
      .iterate((theta, phi)) {
        case (p, t) =>
          leapfrog(t, p, eps)
      }
  }

  private def logAcceptance(
    propTheta: DenseVector[Double],
    propPhi: DenseVector[Double],
    theta: DenseVector[Double],
    phi: DenseVector[Double]) = {

    val ap = ll(propTheta) + propPhi.t * propPhi * 0.5 -
      ll(theta) - phi.t * phi * 0.5
    if (ap.isNaN) {
      -1e99
    } else {
      min(-ap, 0.0)
    }
  }

  /**
    * A single step of the Hamiltonian Monte Carlo Algorithm
    * with dual averaging
    * @param s the current state
    */
  def step(s: HmcState): Rand[HmcState] = {
    for {
      phi <- priorPhi
      (propPhi, propTheta) = leapfrogs(s.theta, phi, eps).
        take(l).last
      a = logAcceptance(propTheta, propPhi, s.theta, phi)
      u <- Uniform(0, 1)
      next = if (log(u) < -a) {
        HmcState(propTheta, s.accepted + 1)
      } else {
        s
      }
    } yield next
  }

  /**
    * Perform HMC for the Gaussian Process Hyper Parameters
    */
  def sample(init: DenseVector[Double]): Process[HmcState] = {
    val initState = HmcState(init, 0)
    MarkovChain(initState)(step)
  }
}

object Hmc {
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
    * Sample GP parameters using HMC with dual averaging
    * @param ys a vector of observations of a GP
    */
  def sampleGp(
    ys:     Vector[Data],
    dist:   (Location[Double], Location[Double]) => Double,
    init:   Parameters,
    ll:     Parameters => Double,
    m:      Int,
    l:      Int,
    eps:    Double) = {

    def newLl = (p: DenseVector[Double]) => ll(vectorToParams(init, p))
    def newGrad =
      (p: DenseVector[Double]) => mllGradient(ys, dist)(vectorToParams(init, p))

    val hmc = Hmc(m, l, eps, newGrad, newLl)

    val theta = paramsToDenseVector(init)
    val initState = HmcState(theta, 0)
    MarkovChain(initState)(hmc.step)
  }
}
