package gp.core

import breeze.stats.distributions._
import breeze.linalg.{DenseVector, diag}
import GaussianProcess._
import math._

/**
  * TODO: How to handle proposals of NaN
  */
case class DualAverageState(
  iter:      Int,
  theta:     DenseVector[Double],
  accepted:  Int,
  logeps:    Double,
  logepsbar: Double,
  hm:        Double)

/**
  * Hamiltonian Monte Carlo with dual averaging for selecting
  * optimum step size
  * @param lambda target simulation length
  * @param delta target acceptance rate
  * @param m the diagonal elements of the multivariate normal
  * prior distribution for the momentum parameters
  * @param mAdapt the number of adaptation iterations to run
  */
case class DualAverage(
  lambda: Double,
  delta: Double,
  m: DenseVector[Double],
  mAdapt: Int,
  gradient: DenseVector[Double] => DenseVector[Double],
  ll: DenseVector[Double] => Double) {

  /**
    * Prior distribution for the momentum variables
    */
  def priorPhi = {
    val zero = DenseVector.zeros[Double](m.size)
    MultivariateGaussian(zero, diag(m))
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
    t0:    Double = 10.0
  )(hm0: Double,
    logeps0: Double,
    logepsbar0: Double): (Double, Double, Double) = {

    val md = m.toDouble
    val ra = 1 / (md + t0)
    val hm = (1 - ra) * hm0 + ra * (delta - acceptProb)
    val logeps1 = mu - (sqrt(md) / gamma) * hm
    val power = pow(md, -k)
    val logepsbar1 = power * logeps1 + (1.0 - power) * logepsbar0

    (hm, logeps1, logepsbar1)
  }

  /**
    * A single step of the Hamiltonian Monte Carlo Algorithm
    * with dual averaging
    * @param s the current state
    */
  def step(mu: Double)(s: DualAverageState): Rand[DualAverageState] = {
    for {
      phi <- priorPhi
      eps = exp(s.logeps)
      lm = max(1, round(lambda / eps).toInt).min(1000)
      (propTheta, propPhi) = Hmc.leapfrogs(eps, gradient, lm)(s.theta, phi)
      a = Hmc.logAcceptance(propTheta, propPhi, s.theta, phi, ll, priorPhi)
      // (hm1, logeps1, logepsbar1) = (s.hm, s.logeps, s.logepsbar)
      (hm1, logeps1, logepsbar1) = if (s.iter < mAdapt) {
        if (s.iter % (mAdapt / 10) == 0)
          println(s"Adaptation Phase: Iteration ${s.iter} / $mAdapt, Epsilon: $eps, Leapfrog steps: $lm")
        updateEps(s.iter, mu, min(1.0, exp(a)))(s.hm, s.logeps, s.logepsbar)
      } else {
        if (s.iter % 1000 == 0)
          println(s"Main sampling run, iteration ${s.iter}, accepted ${s.accepted.toDouble / s.iter}")
        (s.hm, s.logepsbar, s.logepsbar)
      }
      u <- Uniform(0, 1)
      next = if (log(u) < a) {
        // println("accepted")
        // println(s.accepted.toDouble / s.iter)
        DualAverageState(s.iter + 1, propTheta, s.accepted + 1,
          logeps1, logepsbar1, hm1)
      } else {
        s.copy(iter = s.iter + 1, logeps = logeps1,
          logepsbar = logepsbar1, hm = hm1)
      }
    } yield next
  }

  /**
    * Perform HMC for the Gaussian Process Hyper Parameters
    */
  def sample(init: DenseVector[Double]): Process[DualAverageState] = {
    val eps0 = DualAverage.findReasonableEpsilon(init, ll, priorPhi, m, gradient)
    println(s"initial step size $eps0")
    val initState = DualAverageState(1, init, 0, log(eps0), 0.0, 0.0)
    MarkovChain(initState)(step(log(10 * eps0)))
  }
}

object DualAverage {
  /**
    * Initialise the value of epsilon
    */
  def findReasonableEpsilon(
    theta: DenseVector[Double],
    ll: DenseVector[Double] => Double,
    priorPhi: ContinuousDistr[DenseVector[Double]],
    m: DenseVector[Double],
    gradient: DenseVector[Double] => DenseVector[Double]): Double = {
    println("finding reasonable epsilon")

    val eps = 1.0
    val phi = priorPhi.draw
    val (initTheta, initPhi) = Hmc.leapfrog(eps, gradient)(theta, phi)
    def prop(propTheta: DenseVector[Double], propPhi: DenseVector[Double]) =
      Hmc.logAcceptance(propTheta, propPhi, theta, phi, ll, priorPhi)
    val i = prop(initTheta, initPhi) > log(0.5)
    val a = if (i) 1.0 else -1.0

    def loop(thetaP: DenseVector[Double],
      phiP: DenseVector[Double], curEps: Double): Double = {

      if (a * prop(thetaP, phiP) > -a * log(2.0)) {
        val (propTheta, propPhi) = Hmc.leapfrog(curEps, gradient)(theta, phi)
        loop(propTheta, propPhi, pow(2.0, a) * curEps)
      } else {
        curEps
      }
    }

    loop(initTheta, initPhi, eps)
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
    m:      DenseVector[Double],
    lambda: Double,
    delta:  Double,
    mAdapt: Int) = {

    def newLl(p: DenseVector[Double]) =
      ll(Hmc.vectorToParams(init, p))
    def newGrad(p: DenseVector[Double]) =
      mllGradient(ys, dist)(Hmc.vectorToParams(init, p))

    val hmc = DualAverage(lambda, delta, m, mAdapt, newGrad, newLl)
    val priorPhi = MultivariateGaussian(
      DenseVector.zeros[Double](m.size), diag(m))

    val theta = Hmc.paramsToDenseVector(init)
    val eps0 = findReasonableEpsilon(theta, newLl, priorPhi, m, newGrad)
    val initState = DualAverageState(1, theta, 0, log(eps0), 0.0, 0.0)
    MarkovChain(initState)(hmc.step(log(10) + log(eps0)))
  }
}
