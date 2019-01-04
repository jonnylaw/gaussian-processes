package com.github.jonnylaw.gp

import breeze.linalg.{DenseMatrix, DenseVector}
import breeze.stats.distributions._

object Summarise {

  /**
    * Calculate the intervals
    */
  def getInterval(mean: Double, variance: Double, interval: Double) = {
    Gaussian(mean, math.sqrt(variance)).inverseCdf(interval)
  }

  /**
    * Extract the diagonal elements of the Breeze dense vector
    */
  def getDiagonal(m: DenseMatrix[Double]): Vector[Double] = {
    for {
      i <- Vector.range(0, m.cols)
    } yield m(i, i)
  }

  /**
    * Get intervals from marginal Normal distributions
    * @param mean the mean of a multivariate Normal distribution
    * @param covariance the covariance of a multivariate Normal distribution
    * @param interval the desired interval to return between zero and one
    * @return a Vector containing the intervals
    */
  def getIntervals(mean: DenseVector[Double],
                   covariance: DenseMatrix[Double],
                   interval: Double): Vector[(Double, Double)] = {

    (mean.data, getDiagonal(covariance)).zipped.map {
      case (mu, cov) =>
        (getInterval(mu, cov, interval), getInterval(mu, cov, 1.0 - interval))
    }.toVector
  }
}
