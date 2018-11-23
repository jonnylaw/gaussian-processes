package gp.core

import breeze.stats.distributions._
import breeze.linalg.{DenseVector, diag}
import math._

/**
  * No U-turn sampler with dual averaging, HMC with selection of step size
  * epsilon and number of leapfrog steps l in an adaptation phase
  * TODO: This should be refactored to share some code with HMC methods
  * and split into smaller functions
  * @param deltamax if ll(theta) - prior.logPdf(phi) - log(u) < deltamax stop 
  * doubling
  */
case class Nuts(
  m: DenseVector[Double],
  delta: Double,
  mAdapt: Int,
  gradient: DenseVector[Double] => DenseVector[Double],
  ll: DenseVector[Double] => Double,
  deltamax: Double = 1000) {

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
    m:          Int,
    mu:         Double,
    acceptProb: Double,
    nAccept:    Int,
    k:          Double = 0.75,
    gamma:      Double = 0.05,
    t0:         Double = 10
  )(hm0:        Double,
    logeps0:    Double,
    logepsbar0: Double): (Double, Double, Double) = {

    val md = m.toDouble
    val ra = 1 / (md + t0)
    val hm = (1 - ra) * hm0 + ra * (delta - acceptProb / nAccept.toDouble)
    val logeps1 = mu - (sqrt(md) / gamma) * hm
    val power = pow(md, -k)
    val logepsbar1 = power * logeps1 + (1.0 - power) * logepsbar0

    (hm0, logeps1, logepsbar1)
  }

  case class TreeState(
    thetaM: DenseVector[Double],
    phiM: DenseVector[Double],
    thetaP: DenseVector[Double],
    phiP: DenseVector[Double],
    theta1: DenseVector[Double],
    n:      Int,
    s:      Int,
    acceptProb: Double,
    nAccept: Int)

  /**
    * This is very ugly
    */
  def buildTree(
    u:      Double,
    v:      Int,
    j:      Int,
    eps:    Double,
    theta0: DenseVector[Double],
    phi0:   DenseVector[Double])(
    theta:  DenseVector[Double],
    phi:    DenseVector[Double]): TreeState = {

    if (j == 0) {
      val (t1, p1) = Hmc.leapfrog(eps * v, gradient)(theta, phi)
      val a1 = ll(t1) + priorPhi.logPdf(p1)
      val n = if (a1 >= log(u)) 1 else 0
      val s = if (deltamax + a1 > log(u)) 1 else 0
      val a = Hmc.logAcceptance(t1, p1, theta0, phi0, ll, priorPhi)

      TreeState(t1, p1, t1, p1, t1, n, s, a, 1)
    } else {
      val bt = buildTree(u, v, j - 1, eps, theta0, phi0) _
      val st1 = bt(theta, phi)
      val st = if (st1.s == 1) {
        val st2 = if (v == -1) {
          bt(st1.thetaM, st1.phiM)
        } else {
          bt(st1.thetaP, st1.phiP)
        }
        val u = Uniform(0, 1).draw
        val p = st2.n.toDouble / (st1.n.toDouble + st2.n.toDouble)
        val newTheta = if (u < p) {
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

    // println("updating S")
    val i1 = if (((thetaP - thetaM) dot phiM) >= 0) 1 else 0
    val i2 = if (((thetaP - thetaM) dot phiP) >= 0) 1 else 0
    val res = s * i1 * i2
    // println(s"s is $res")
    res
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
        // println("Current s is one")
        val newState = for {
          vj <- sampleDirection
          // _ = println(s"direction $vj")
          bt = buildTree(u, vj, j, eps, theta, phi0) _
          st1 = if (vj == -1) {
            bt(current.thetaM, current.phiM)
          } else {
            bt(current.thetaP, current.phiP)
          }
          u <- Uniform(0, 1)
          p = st1.n.toDouble / current.n.toDouble
          // _ = println(s"s prime is ${st1.s}")
          // _ = println(s"acceptance prob n prime / n is $p")
          theta = if (st1.s == 1 && u < p) {
            // println("accepted!")
            st1.theta1
          } else {
            current.theta1
          }
        } yield st1.copy(theta1 = theta, n = current.n + st1.n,
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
  def step(mu: Double)(s: DualAverageState): Rand[DualAverageState] = {
    for {
      phi <- priorPhi
      u <- Uniform(0, exp(ll(s.theta) + priorPhi.logPdf(phi)))
      eps = exp(s.logeps)
      initst = TreeState(s.theta, phi, s.theta, phi, s.theta, 1, 1, 0.0, 0)
      st = loopTrees(u, eps, 0, phi, s.theta)(initst)
      // (hm1, logeps1, logepsbar1) = (s.hm, s.logeps, s.logepsbar)
      (hm1, logeps1, logepsbar1) = if (s.iter < mAdapt) {
        if (s.iter % (mAdapt / 10) == 0)
          println(s"Optimising step size and no. leapfrogs, iteration ${s.iter} / $mAdapt, current step size $eps")
        updateEps(s.iter, mu, min(1.0, exp(st.acceptProb)), st.nAccept)(s.hm, s.logeps, s.logepsbar)
      } else {
        if (s.iter % 1000 == 0)
          println(s"Main sampling run, iteration ${s.iter}")
        (s.hm, s.logepsbar, s.logepsbar)
      }
      next = DualAverageState(s.iter + 1, st.theta1,
        s.accepted + st.nAccept, logeps1, logepsbar1, hm1)
    } yield next
  }

  /**
    * Perform HMC for the Gaussian Process Hyper Parameters
    */
  def sample(init: DenseVector[Double]): Process[DualAverageState] = {
    val eps0 = DualAverage.
      findReasonableEpsilon(init, ll, priorPhi, m, gradient)
    // println(s"initial step size $eps0")
    val initState = DualAverageState(1, init, 0, log(eps0), 0.0, 0.0)
    MarkovChain(initState)(step(log(10 * eps0)))
  }
}
