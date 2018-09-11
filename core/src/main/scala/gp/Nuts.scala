package gp.core

import breeze.stats.distributions._
import breeze.linalg.{DenseVector, diag}
import math._

/**
  * No U-turn sampler with dual averaging, HMC with selection of step size
  * epsilon and number of leapgrog steps l in an adaptation phase
  * TODO: This should be refactored to share some code with HMC methods
  * and split into smaller functions
  * TODO: What is deltamax!?
  */
case class Nuts(
  delta: Double,
  m: DenseVector[Double],
  mAdapt: Int,
  gradient: DenseVector[Double] => DenseVector[Double],
  ll: DenseVector[Double] => Double,
deltamax: Double) {

  def findReasonableEpsilon(theta: DenseVector[Double]): Double = {
    val eps = 1.0
    val phi = priorPhi.draw
    val (initTheta, initPhi) = leapfrog(phi, theta, eps)
    val prop = (propTheta: DenseVector[Double], propPhi: DenseVector[Double]) =>
      logAcceptance(propTheta, propPhi, phi, theta)
    val i = prop(initTheta, initPhi) > log(0.5)
    val a = if (i) 2 else -1

    def loop(thetaP: DenseVector[Double],
      phiP: DenseVector[Double], curEps: Double): Double = {

      if (pow(prop(thetaP, phiP), a) > pow(2, -a)) {
        val (propTheta, propPhi) = leapfrog(phiP, thetaP, eps)
        loop(propTheta, propPhi, pow(2, a) * curEps)
      } else {
        curEps
      }
    }

    loop(initTheta, initPhi, eps)
  }

  private def leapfrogStep(
    phi: DenseVector[Double],
    theta: DenseVector[Double],
    eps: Double) = {

    val newTheta = theta + eps * diag(m.map(1.0 / _)) * phi
    val newPhi = phi + eps * gradient(newTheta)

    (newPhi, newTheta)
  }

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

  private def logAcceptance(
    propTheta: DenseVector[Double],
    propPhi: DenseVector[Double],
    theta: DenseVector[Double],
    phi: DenseVector[Double]) = {

    val ap = ll(propTheta) + priorPhi.logPdf(propPhi) -
      ll(theta) - priorPhi.logPdf(phi)
    if (ap.isNaN) {
      -1e99
    } else {
      (-ap).min(0.0)
    }
  }

  /**
    * Perform a leapfrog update with 1 step
    */
  def leapfrog(
    theta: DenseVector[Double],
    phi: DenseVector[Double],
    eps: Double) = {

    val (p1) = leapfrogHalfStep(theta, phi, eps)
    val (t1, p2) = leapfrogStep(theta, p1, eps)
    val p3 = leapfrogHalfStep(t1, p2, eps)

    (p3, t1)
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

    val hm = (1 - 1 / (m + t0)) * hm0 + 1 / (m + t0) * (delta - acceptProb)
    val logem = mu - (sqrt(m) / gamma) * hm
    val logem1 = pow(m, -k) * logem + (1 - pow(m, -k)) * logeps0

    (hm0, logem1)
  }

  case class TreeState(
    thetaP: DenseVector[Double],
    phiP: DenseVector[Double],
    thetaM: DenseVector[Double],
    phiM: DenseVector[Double],
    theta1: DenseVector[Double],
    n: Int,
    s: Int,
    acceptProb: Double,
    nAccept: Int)
    
  /**
    * This is very ugly
    */
  def buildTree(
    u:        Double,
    v:        Int,
    j:        Int,
    eps:      Double,
    theta0: DenseVector[Double],
    phi0: DenseVector[Double])(
    theta: DenseVector[Double],
    phi: DenseVector[Double]): TreeState = { 

    if (j == 0) {
      val (p1, t1) = leapfrog(theta, phi, v * eps)
      val n = if (ll(t1) + priorPhi.logPdf(p1) > log(u)) 1 else 0
      val s = if (deltamax + ll(t1) + priorPhi.logPdf(p1) > log(u)) 1 else 0
      val a = logAcceptance(t1, p1, theta0, phi0)
      TreeState(t1, p1, t1, p1, t1, n, s, min(1, exp(a)), 1)
    } else {
      // recursively build subtrees
      val bt = buildTree(u, v, j - 1, eps, theta0, phi0) _
      val st1 = bt(theta, phi)
      val st = if (st1.s == 1) {
        val st2 = if (v == -1) {
          bt(st1.thetaM, st1.phiM)
        } else {
          bt(st1.thetaP, st1.phiP)
        }
        val u = scala.util.Random.nextDouble()
        val newTheta = if (u < st2.n / (st1.n + st2.n)) {
          st2.theta1
        } else {
          st1.theta1
        }
        val a = st1.acceptProb + st2.acceptProb
        val nA = st1.nAccept + st2.nAccept
     val newS = updateS(st2.s, st2.thetaP, st2.thetaM, st2.phiP, st2.phiM)
        st2.copy(theta1 = newTheta, acceptProb = a, nAccept = nA,
          s = newS, n = st1.n + st2.n)
      } else {
        st1
      }

      st
    }
  }

  def updateS(
    s: Int,
    thetaP: DenseVector[Double],
    thetaM: DenseVector[Double],
    phiP: DenseVector[Double],
    phiM: DenseVector[Double]) = {

    val i1 = if (((thetaP - thetaM) dot phiM) >= 0) 1 else 0
    val i2 = if (((thetaP - thetaM) dot phiP) >= 0) 1 else 0
    s * i1 * i2
  }

  def sampleDirection: Rand[Int] = {
    val b = new Bernoulli(0.5) 

    b.map(e => if (e) -1 else 1)
  }

  /**
    * Urgh!
    */
  def loopTrees(
    u: Double,
    eps: Double,
    j: Int,
    phi0: DenseVector[Double],
    theta: DenseVector[Double]) = { st: TreeState =>

    def loop(current: TreeState, j: Int): TreeState = {
      if (current.s == 1) {
        val newState = for {
          vj <- sampleDirection
          bt = buildTree(u, vj, j, eps, theta, phi0) _
          st1 = if (vj == -1) {
            bt(st.thetaM, st.phiM)
          } else {
            bt(st.thetaP, st.phiP)
          }
          u <- Uniform(0, 1)
          theta = if (st1.s == -1 && u < st1.n / st.n) {
            st1.theta1
          } else {
            st.theta1
          }
        } yield st1.copy(theta1 = theta, n = st.n + st1.n,
            s = updateS(st1.s, st1.thetaP, st1.thetaM, st1.phiP, st1.phiM))
        loop(newState.draw, j + 1)
      } else {
        current
      }
    }

    loop(st, j)
  }

  /**
    * A single step of the Hamiltonian Monte Carlo Algorithm
    * @param s the current state
    */
  def step(s: HmcState): Rand[HmcState] = {
    for {
      phi <- priorPhi
      u <- Uniform(0, exp(ll(s.theta) + priorPhi.logPdf(phi)))
      eps = exp(s.logeps)
      initst = TreeState(s.theta, phi, s.theta, phi, s.theta, 1, 1, 0.0, 0)
      st = loopTrees(u, eps, 0, phi, s.theta)(initst)
      (hm1, logeps1) = if (s.iter < mAdapt) {
        updateEps(s.iter, log(10) + s.logeps,
          min(1.0, exp(st.acceptProb)))(s.hm, s.logeps)
      } else {
        (s.hm, s.logeps)
      }
      next = HmcState(s.iter + 1, st.theta1,
        s.accepted + st.nAccept, logeps1, hm1)
    } yield next
  }
}
