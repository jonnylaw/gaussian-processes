package com.github.jonnylaw.gp

import scala.annotation.tailrec
import breeze.linalg.{DenseMatrix, DenseVector}
import breeze.stats.distributions.{Rand, MultivariateGaussian}

object Hmc {
  def priorPhi(m: DenseMatrix[Double]) = {
    val zero = DenseVector.zeros[Double](m.cols)
    MultivariateGaussian(zero, m).map(_.data)
  }

  def phiHalfStep(eps: Double,
                  theta: Array[Double],
                  phi: Array[Double],
                  gradient: Array[Double] => Array[Double]) =
    phi.zip(gradient(theta)).map { case (p, dt) => p + dt * 0.5 * eps }

  def thetaStep(eps: Double, theta: Array[Double], phi: Array[Double]) =
    theta.zip(phi).map { case (t, p) => t + eps * p }

  /**
    * Perform a leapfrog full step with the momentum
    * @param theta the current value of the parameters
    * @param phi the current value of the momentum
    * @param eps the step size
    * @return a tuple containing the updated value of the parameters and
    * momentum
    */
  def leapfrog(eps: Double,
               gradient: Array[Double] => Array[Double],
               theta: Array[Double],
               phi: Array[Double]) = {
    val p1 = phiHalfStep(eps, theta, phi, gradient)
    val t1 = thetaStep(eps, theta, p1)
    val p2 = phiHalfStep(eps, t1, p1, gradient)

    (t1, p2)
  }

  /**
    * Perform l leapfrog steps
    */
  @tailrec
  def leapfrogs(eps: Double,
                gradient: Array[Double] => Array[Double],
                l: Int,
                theta: Array[Double],
                phi: Array[Double]): (Array[Double], Array[Double]) = {
    if (l == 0) {
      (theta, phi)
    } else {
      val (t, p) = leapfrog(eps, gradient, theta, phi)
      leapfrogs(eps, gradient, l - 1, t, p)
    }
  }

  // E_k = 0.5 p^2
  def kinetic(phi: Array[Double]): Double = {
    var k = 0.0
    var i = 0
    while (i < phi.size) {
      val p = phi(i)
      k += (p * p)
      i += 1
    }
    k * 0.5
  }

  def logDensity(theta: Array[Double],
                 phi: Array[Double],
                 pos: Array[Double] => Double) =
    pos(theta) - kinetic(phi)

  /**
    * Calculate the log-acceptance rate
    */
  def logAcceptance(propTheta: Array[Double],
                    propPhi: Array[Double],
                    theta: Array[Double],
                    phi: Array[Double],
                    pos: Array[Double] => Double) = {
    val a = logDensity(propTheta, propPhi, pos) - logDensity(theta, phi, pos)
    if (a.isNaN) { Math.log(0.0) } else { a.min(0.0) }
  }

  def step(eps: Double,
           l: Int,
           m: DenseMatrix[Double],
           pos: Array[Double] => Double,
           gradient: Array[Double] => Array[Double])(
    theta: Array[Double]): Rand[(Array[Double], Double)] = {

    for {
      phi <- priorPhi(m)
      (t, p) = Hmc.leapfrogs(eps, gradient, l, theta, phi)
      a = Hmc.logAcceptance(t, p, theta, phi, pos)
      u <- breeze.stats.distributions.Uniform(0, 1)
      next = if (math.log(u) < a) {
        (t, a)
      } else {
        (theta, a)
      }
    } yield next
  }
}
