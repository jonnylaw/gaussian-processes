package gp.core

import breeze.stats.distributions._
import breeze.linalg.{DenseVector, diag}
import GaussianProcess._
import math._

/**
  * TODO: How to handle proposals of NaN
  */
case class HmcDaState(
  iter:     Int,
  theta:    DenseVector[Double],
  accepted: Int,
  logeps:   Double,
  hm:       Double)

/**
  * Hamiltoniam Monte Carlo with dual averaging for selecting
  * optimum step size
  * @param lambda target simulation length
  * @param delta target acceptance rate
  * @param m the diagonal elements of the multivariate normal
  * prior distribution for the momentum parameters
  * @param mAdapt the number of adaptation iterations to run
  */
case class HmcDa(
  lambda: Double,
  delta: Double,
  m: DenseVector[Double],
  mAdapt: Int,
  gradient: DenseVector[Double] => DenseVector[Double],
  ll: DenseVector[Double] => Double) {

  private def leapfrogHalfStep(
    theta: DenseVector[Double],
    phi: DenseVector[Double],
    eps: Double) = 
      phi + eps * 0.5 * gradient(theta)

  /**
    * Prior distribution for the momentum variables
    */
  def priorPhi = {
    val zero = DenseVector.zeros[Double](m.size)
    MultivariateGaussian(zero, diag(m))
  }

  /**
    * Perform a leapfrog update with 1 step
    */
  def leapfrog(
    theta: DenseVector[Double],
    phi: DenseVector[Double],
    eps: Double) = {

    val p1 = leapfrogHalfStep(theta, phi, eps)
    val t1 = theta + eps * diag(m.map(1.0 / _)) * p1
    val p3 = leapfrogHalfStep(t1, p1, eps)

    (t1, p3)
  }

  private def leapfrogs(
    theta: DenseVector[Double],
    phi: DenseVector[Double],
    eps: Double): Stream[(DenseVector[Double], DenseVector[Double])] = {

    Stream
      .iterate((theta, phi)) {
        case (t, p) =>
          leapfrog(t, p, eps)
      }
  }

  private def logAcceptance(
    propTheta: DenseVector[Double],
    propPhi: DenseVector[Double],
    theta: DenseVector[Double],
    phi: DenseVector[Double]) = {

    val ap = ll(propTheta) + priorPhi.logPdf(propPhi) +
      ll(theta) - priorPhi.logPdf(phi)

    if (ap.isNaN) {
      -1e99
    } else {
      ap
    }
  }

  def findReasonableEpsilon(theta: DenseVector[Double]): Double = {
    println("finding reasonable epsilon")

    val eps = 1.0
    val phi = priorPhi.draw
    val (initTheta, initPhi) = leapfrog(theta, phi, eps)
    def prop(propTheta: DenseVector[Double], propPhi: DenseVector[Double]) =
      logAcceptance(propTheta, propPhi, theta, phi)
    val i = prop(initTheta, initPhi) > log(0.5)
    val a = if (i) 1.0 else -1.0

    def loop(thetaP: DenseVector[Double],
      phiP: DenseVector[Double], curEps: Double): Double = {

      if (a * prop(thetaP, phiP) > -a * log(2)) {
        val (propTheta, propPhi) = leapfrog(theta, phi, curEps)
        loop(propTheta, propPhi, pow(2, a) * curEps)
      } else {
        curEps
      }
    }

    loop(initTheta, initPhi, eps)
  }

  /**
    * Update the value of epsilon during an adaptation phase
    * @param m the iteration number
    * @param mu default set to log(10 * eps0)
    * @param acceptProb the acceptance probability from the previous time step
    */
  def updateEps(
    m:     Int,
    mu:    Double,
    acceptProb: Double,
    k:     Double = 0.75,
    gamma: Double = 0.05,
    t0:    Double = 10
    )(hm0: Double, logeps0: Double): (Double, Double) = {

    val ra = 1 / (m + t0)
    val hm = (1 - ra) * hm0 + ra * (delta - acceptProb)
    val logem = mu - (sqrt(m) / gamma) * hm
    val logem1 = pow(m, -k) * logem + (1 - pow(m, -k)) * logeps0

    (hm, logem1)
  }

  /**
    * A single step of the Hamiltonian Monte Carlo Algorithm
    * with dual averaging
    * @param s the current state
    */
  def step(s: HmcDaState): Rand[HmcDaState] = {
    for {
      phi <- priorPhi
      eps = exp(s.logeps) // use current value of step size
      lm = max(1, round(lambda / eps).toInt)
      (propPhi, propTheta) = leapfrogs(s.theta, phi, eps).
        take(lm).last
      _ = println(s"proposed parameters $propTheta")
      a = logAcceptance(propTheta, propPhi, s.theta, phi)
      (hm1, logeps1) = if (s.iter < mAdapt) {
        updateEps(s.iter, log(10) + s.logeps,
          min(1.0, exp(a)))(s.hm, s.logeps)
      } else {
        (s.hm, s.logeps)
      }
      _ = if (s.iter < 1000 && s.iter % 100 == 0) {
        println(s"current iteration: ${s.iter}")
        println(s"current step size: ${exp(logeps1)}")
        println(s"current hm: ${hm1}")
        println(s"current number of steps $lm")
      }
      u <- Uniform(0.0, 1.0)
      next = if (log(u) < a) {
        HmcDaState(s.iter + 1, propTheta, s.accepted + 1, logeps1, hm1)
      } else {
        s.copy(iter = s.iter + 1, logeps = logeps1, hm = hm1)
      }
    } yield next
  }

  /**
    * Perform HMC for the Gaussian Process Hyper Parameters
    */
  def sample(init: DenseVector[Double]): Process[HmcDaState] = {
    val eps0 = findReasonableEpsilon(init)
    println(s"initial step size $eps0")
    val initState = HmcDaState(1, init, 0, log(eps0), 0.0)
    MarkovChain(initState)(step)
  }
}

object HmcDa {
  /**
    * Sample GP parameters using HMC with dual averaging
    * @param ys a vector of observations of a GP
    */
  def sampleGp(
    ys:     Vector[Data],
    dist:   (Location[Double], Location[Double]) => Double,
    init:   Parameters,
    ll:     Parameters => Double,
    m:      DenseVector[Double],
    lambda: Double,
    delta:  Double,
    mAdapt: Int) = {

    def newLl = (p: DenseVector[Double]) => ll(Hmc.vectorToParams(init, p))
    def newGrad =
      (p: DenseVector[Double]) => mllGradient(ys, dist)(Hmc.vectorToParams(init, p))

    val hmc = HmcDa(lambda, delta, m, mAdapt, newGrad, newLl)

    val theta = Hmc.paramsToDenseVector(init)
    val eps0 = hmc.findReasonableEpsilon(theta)
    val initState = HmcDaState(1, theta, 0, log(eps0), 0.0)
    MarkovChain(initState)(hmc.step)
  }
}
