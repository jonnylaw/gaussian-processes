package gp.core

import breeze.stats.distributions._
import breeze.linalg.{DenseVector, DenseMatrix}
import breeze.numerics.log1p
import GaussianProcess._
import math._

/**
  * Discrete uniform distribution
  */
case class DiscreteUniform(min: Int, max: Int) extends Rand[Int] {
  def draw =
    min + scala.util.Random.nextInt(max - min + 1)
}

case class HmcState(
  iter:     Int,
  theta:    DenseVector[Double],
  accepted: Int)

/**
  * Hamiltonian Monte Carlo on unconstrained parameter spaces
  * @param m the dimension of the parameter space
  * @param l the number of leapgfrog steps
  * @param eps the step size
  * @param gradient
  * @param ll 
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
    * A single step of the Hamiltonian Monte Carlo Algorithm
    * with dual averaging
    * @param s the current state
    */
  def step(s: HmcState): Rand[HmcState] = {
    for {
      phi <- priorPhi
      // _ = println(s"initial momentum $phi")
      (propTheta, propPhi) = Hmc.leapfrogs(eps, gradient, l)(s.theta, phi)
      // _ = println(s"proposed theta after $l leapfrogs: $propTheta")
      // _ = println(s"proposed phi: $propPhi")
      a = Hmc.logAcceptance(propTheta, propPhi, s.theta, phi, ll, priorPhi)
      u <- Uniform(0, 1)
      _ = if (s.iter % 1000 == 0) {
        println(s"Main sampling run, iteration ${s.iter}, accepted ${s.accepted.toDouble / s.iter}")
      }
      // _ = println(s"u is ${log(u)}")
      // _ = println(s"a is ${a}")
      next = if (log(u) < a) {
        // println("accepted!")
        HmcState(s.iter + 1, propTheta, s.accepted + 1)
      } else {
        s.copy(iter = s.iter + 1)
      }
    } yield next
  }

  /**
    * Perform HMC for the Gaussian Process Hyper Parameters
    */
  def sample(init: DenseVector[Double]): Process[HmcState] = {
    val initState = HmcState(0, init, 0)
    MarkovChain(initState)(step)
  }
}

object Hmc {
  /**
    * Perform a leapfrog full step with the momentum
    * @param theta the current value of the parameters
    * @param phi the current value of the momentum
    * @param eps the step size
    * @return a tuple containing the updated value of the parameters and
    * momentum 
    */
  def leapfrog(
    eps: Double,
    gradient: DenseVector[Double] => DenseVector[Double])(
    theta: DenseVector[Double],
    phi: DenseVector[Double]) = {
    val p1 = phi - eps * 0.5 * gradient(theta)
    // println(s"momentum is $p1")
    val t1 = theta + eps * p1
    // println(s"parameters are $t1, constrained ${constrain(t1).map(_.constrained)}")
    val p2 = p1 - eps * 0.5 * gradient(t1)

    (t1, p2)
  }

  /**
    * Perform l leapfrog steps
    */
  def leapfrogs(
    eps: Double,
    gradient: DenseVector[Double] => DenseVector[Double],
    l: Int)(
    theta: DenseVector[Double],
    phi: DenseVector[Double]): (DenseVector[Double], DenseVector[Double]) = 
    (1 to l).foldLeft((theta, phi)){ case ((t, p), _) =>
      leapfrog(eps, gradient)(t, p)
    }

  def constrain(p: DenseVector[Double]): Vector[Double] = {
    val phi = Hmc.bounded(0.0, 1.0)(p(0))
    val mu = Hmc.unbounded(p(1))
    val sigma = Hmc.boundedBelow(0)(p(2))

    Vector(phi.constrained, mu.constrained, sigma.constrained)
  }

  /**
    * Calculate the log-acceptance rate
    */
  def logAcceptance(
    propTheta: DenseVector[Double],
    propPhi: DenseVector[Double],
    theta: DenseVector[Double],
    phi: DenseVector[Double],
    ll: DenseVector[Double] => Double,
    priorPhi: ContinuousDistr[DenseVector[Double]]) = {

    // check if any parameters are NaN
    val thetaNans = propTheta.data.map(_.isNaN).
      foldLeft(false)((f, a) => a || f)

    // println(s"prop theta ${constrain(propTheta)}")

    if (thetaNans) {
      -1e99
    } else {
      val ap = ll(propTheta) + priorPhi.logPdf(propPhi) -
      ll(theta) - priorPhi.logPdf(phi)
      // println(s"proposed momentum: $propPhi")
      // println(s"Log likelihood of theta ${constrain(theta)} is ${ll(theta)}")
      // println(s"Log likelihood of prop theta ${constrain(propTheta)} is ${ll(propTheta)}")
      // println(s"the acceptance probability is $ap")

      if (ap.isNaN) {
        -1e99
      } else {
        ap
      }
    }
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
    * Sample GP parameters using HMC with dual averaging
    * @param ys a vector of observations of a GP
    */
  def sampleGp(
    ys:   Vector[Data],
    dist: (Location[Double], Location[Double]) => Double,
    init: Parameters,
    ll:   Parameters => Double,
    m:    Int,
    l:    Int,
    eps:  Double
  ) = {

    def newLl(p: DenseVector[Double]) = {
      val params = vectorToParams(init, p)
      val kp = KernelParameters.constrainParams(params.kernelParameters)
      ll(params.copy(kernelParameters = kp))
    }
    def newGrad(p: DenseVector[Double]) = {
      val params = vectorToParams(init, p)
      val kp = KernelParameters.constrainParams(params.kernelParameters)
      mllGradient(ys, dist)(params.copy(kernelParameters = kp))
    }

    val kp = KernelParameters.unconstrainParams(init.kernelParameters)
    val theta = paramsToDenseVector(init.copy(kernelParameters = kp))
    val initState = HmcState(0, theta, 0)
    MarkovChain(initState)(Hmc(m, l, eps, newGrad, newLl).step)
  }

  /**
    * Perform a logistic transformation 
    */
  def logistic(x: Double): Double =
    1.0 / (1.0 + exp(-x))

  def logit(p: Double): Double = 
    log(p / (1 - p))

  def softplus(x: Double): Double =
    log1p(exp(x))

  // def boundedBelow1(min: Double)(x: Double) = Parameter(
  //   softplus _,
  //   x => log(exp(x) - 1),
  //   x => x + log(exp(x) - 1),
  //   x => 1 - exp(x) / (exp(x) - 1)
  // )

  /**
    * Parameter in an HMC algorithm
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
