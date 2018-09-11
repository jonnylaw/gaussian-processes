package gp.core

import breeze.stats.distributions._
import breeze.linalg.{DenseMatrix, DenseVector}
import GaussianProcess._
import cats.implicits._

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

case class Matern(sigma: Double, nu: Double, l: Double)
    extends KernelParameters {
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
    * Sample kernel hyperparameters using the metropolis-hastings algorithm
    * @param obs the observations
    * @param dist a distance function
    * @param prior the prior distribution of the hyperparameters
    * @param prop a symmetric proposal distribution
    */
  def sample(
      obs: Vector[Data],
      dist: (Location[Double], Location[Double]) => Double,
      prior: Vector[KernelParameters] => Double,
      prop: Vector[KernelParameters] => Rand[Vector[KernelParameters]]) = {

    val proposal = (p: Parameters) =>
      for {
        kp <- prop(p.kernelParameters)
      } yield p.copy(kernelParameters = kp)

    val ll = (p: Parameters) =>
      prior(p.kernelParameters) +
        GaussianProcess.loglikelihood(obs, dist)(p)

    MarkovChain.Kernels.metropolis(proposal)(ll)
  }

  /**
    * Calculate the gradient of the Kernel Parameters for a given distance
    */
  def gradient(p: KernelParameters)(dist: Double): Vector[Double] = p match {
    case SquaredExp(h, s) =>
      Vector(math.exp(-(dist * dist) / (s * s)),
             -(2 * h / (s * s * s)) * math.exp(-(dist * dist) / (s * s)))
    case White(s)        => Vector(1.0)
    case Matern(_, _, _) => throw new Exception("Not implemented yet")
  }

  /**
    * Build the tangent matrices (a matrix of partial derivatives)
    */
  def tangentMatrix(
      observed: Vector[Data],
      dist: (Location[Double], Location[Double]) => Double,
      p: Vector[KernelParameters]): Vector[DenseMatrix[Double]] = {

    val xs = observed.map(_.x)
    val m = distanceMatrix(xs, dist)
    val res = m.data.toVector map (d => p.flatMap(ps => gradient(ps)(d)))

    res.transpose.map(ms => new DenseMatrix(m.rows, m.cols, ms.toArray))
  }

  // TODO: This is probably a fold, dropping elements from pActual to the next // element of the vector p
  def vectorToParams(p: Vector[KernelParameters],
                     pActual: Vector[Double]): Vector[KernelParameters] = {

    p.foldLeft((Vector[KernelParameters](), pActual)) {
        case ((acc, pa), ps) =>
          ps match {
            case _: SquaredExp =>
              (acc :+ se(pa.head, pa(1)), pa.drop(2))
            case _: White =>
              (acc :+ white(pa.head), pa.tail)
            case _: Matern =>
              (acc :+ matern(pa.head, pa(1), pa(2)), pa.drop(3))
          }
      }
      ._1

  }
}
