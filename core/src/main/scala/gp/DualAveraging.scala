package com.github.jonnylaw.gp

import breeze.stats.distributions._
import breeze.linalg.DenseMatrix
import math._

case class DualAverageState(iter: Int,
                            theta: Array[Double],
                            logeps: Double,
                            logepsbar: Double,
                            hm: Double)

/**
  * Hamiltonian Monte Carlo with dual averaging for selecting
  * optimum step size
  * @param lambda target simulation length
  * @param delta target acceptance rate
  */
object DualAverage {
  /**
    * Update the value of epsilon during an adaptation phase
    * @param m the iteration number
    * @param mu default set to log(10 * eps0)
    * @param acceptProb the acceptance probability from the previous time step
    */
  def updateEps(
      m: Int,
      mu: Double,
      delta: Double,
      acceptProb: Double,
      k: Double = 0.75,
      gamma: Double = 0.05,
      t0: Double = 10.0
  )(hm0: Double,
    logeps0: Double,
    logepsbar0: Double): (Double, Double, Double) = {

    val md = m.toDouble
    val ra = 1 / (md + t0)
    val hm = (1 - ra) * hm0 + ra * (delta - acceptProb)
    val logeps1 = mu - ((sqrt(md) * hm) / gamma)
    val power = pow(md, -k)
    val logepsbar1 = power * logeps1 + (1.0 - power) * logepsbar0

    (hm, logeps1, logepsbar1)
  }

  def step(l0: Int,
           mu: Double,
           m: DenseMatrix[Double],
           pos: Array[Double] => Double,
           gradient: Array[Double] => Array[Double],
           delta: Double = 0.65)(s: DualAverageState) = {

    val eps = exp(s.logeps)
    for {
      (t, a) <- Hmc.step(eps, l0, m, pos, gradient)(s.theta)
      acceptProb = min(1.0, exp(a))
      (hm, logEps, logEpsBar) =
      DualAverage.updateEps(s.iter, mu, delta, acceptProb)(s.hm,
                                                           s.logeps,
                                                           s.logepsbar)

    } yield DualAverageState(s.iter + 1, t, hm, logEps, logEpsBar)
  }

  /**
    * Initialise the value of epsilon
    */
  def findReasonableEpsilon(
      theta: Array[Double],
      phi: Array[Double],
      pos: Array[Double] => Double,
      gradient: Array[Double] => Array[Double]): Double = {

    val eps = 1.0
    val (initTheta, initPhi) = Hmc.leapfrog(eps, gradient, theta, phi)
    def prop(propTheta: Array[Double], propPhi: Array[Double]) =
      Hmc.logAcceptance(propTheta, propPhi, theta, phi, pos)
    val i = prop(initTheta, initPhi) > log(0.5)
    val a = if (i) 1.0 else -1.0

    def loop(thetaP: Array[Double],
             phiP: Array[Double],
             curEps: Double,
             count: Int): Double = {

      if (a * prop(thetaP, phiP) > -a * log(2.0) && count < 100) {
        val (propTheta, propPhi) = Hmc.leapfrog(curEps, gradient, theta, phi)
        loop(propTheta, propPhi, pow(2.0, a) * curEps, count + 1)
      } else if (count > 100) {
        println("Could not find reasonable epsilon in 100 steps")
        curEps
      } else {
        curEps
      }
    }

    loop(initTheta, initPhi, eps, 0)
  }

  def tuneStepsize(k: Int,
                   params: Array[Double],
                   l0: Int,
                   mu: Double,
                   m: DenseMatrix[Double],
                   pos: Array[Double] => Double,
                   gradient: Array[Double] => Array[Double],
                   delta: Double = 0.65): Double = {

    val phi = Hmc.priorPhi(m).draw
    val eps0 = DualAverage.findReasonableEpsilon(params, phi, pos, gradient)
    val init = DualAverageState(1, params, log(eps0), 0.0, 0.0)

    MarkovChain(init)(step(l0, mu, m, pos, gradient)).
      steps.
      drop(k).
      next.
      logepsbar
  }
}
