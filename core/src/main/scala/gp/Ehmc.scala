package com.github.jonnylaw.gp

import cats.implicits._
import math._
import breeze.stats.distributions._
import breeze.linalg.DenseMatrix

object Ehmc {
  def isUTurn(
    theta: Array[Double],
    thetaPlus: Array[Double],
    phiPlus: Array[Double]): Boolean = {

    val nVars = theta.size
    var out = 0.0
    var i = 0
    while (i < nVars) {
      out += (thetaPlus(i) - theta(i)) * phiPlus(i)
      i += 1
    }

    if (out.isNaN)
      true
    else
      out < 0
  }

  def longestBatch(theta: Array[Double],
                   phi: Array[Double],
                   eps: Double,
                   gradient: Array[Double] => Array[Double],
                   L: Int) = {

    def loop(thetaOut: Array[Double],
             phiOut: Array[Double],
             t: Array[Double],
             p: Array[Double],
             l: Int): (Array[Double], Array[Double], Int) = {
      if (!isUTurn(t, p, theta)) {
        val (t1, p1) = Hmc.leapfrog(eps, gradient, t, p)
        if (l == L)
          loop(t1, p1, t1, p1, l + 1)
        else
          loop(thetaOut, phiOut, t1, p1, l + 1)
      } else {
        (thetaOut, phiOut, l)
      }
    }

    loop(theta, phi, theta, phi, 1)
  }

  def longestBatchStep(
      eps: Double,
      l0: Int,
      m: DenseMatrix[Double],
      gradient: Array[Double] => Array[Double],
      pos: Array[Double] => Double
  )(thetal: (Array[Double], Int)): Rand[(Array[Double], Int)] = {

    for {
      phi <- Hmc.priorPhi(m)
      theta = thetal._1
      (t, p, l) = longestBatch(theta, phi, eps, gradient, l0)
      (propTheta, propPhi) = if (l < l0) {
        Hmc.leapfrogs(eps, gradient, l0 - l, t, p)
      } else {
        (t, p)
      }
      u <- Uniform(0, 1)
      a = Hmc.logAcceptance(propTheta, propPhi, theta, phi, pos)
      next = if (log(u) < a) {
        (propTheta, l)
      } else {
        (theta, l)
      }
    } yield next
  }

  def discreteUniform(min: Int, max: Int) =
    for {
      u <- Uniform(0, 1)
    } yield math.floor(u * (max - min) + min).toInt

  def empiricalLongestStep(
    eps: Double,
    l0: Int,
    m: DenseMatrix[Double],
    pos: Array[Double] => Double,
    gradient: Array[Double] => Array[Double],
    theta: Array[Double],
    k: Int = 2000): Vector[Int] = {

    MarkovChain((theta, l0))(longestBatchStep(eps, l0, m, gradient, pos))
      .steps
      .take(k)
      .toVector
      .map(_._2)
  }

  def sample(
    l0: Int,
    mu: Double,
    m: DenseMatrix[Double],
    pos: Array[Double] => Double,
    gradient: Array[Double] => Array[Double],
    warmupIterations: Int,
    initTheta: Array[Double],
    delta: Double = 0.65,
    k: Int = 2000): Process[Array[Double]] = {

    val m = DenseMatrix.eye[Double](initTheta.size)
    val eps = DualAverage.tuneStepsize(warmupIterations, initTheta, l0,
                                       mu, m, pos, gradient, delta)
    val empiricalL = empiricalLongestStep(eps, l0, m, pos, gradient, initTheta, k)

    MarkovChain(initTheta) { theta =>
      for {
        i <- discreteUniform(0, k)
        l = empiricalL(i)
        (newTheta, _) <- Hmc.step(eps, l, m, pos, gradient)(theta)
      } yield newTheta }
  }
}
