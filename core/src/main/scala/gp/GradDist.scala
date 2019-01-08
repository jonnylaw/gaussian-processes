package com.github.jonnylaw.gp

import breeze.stats.distributions._

abstract class GradDist[A](dist: ContinuousDistr[A]) extends ContinuousDistr[A] {
  def logNormalizer: Double =
    dist.logNormalizer
  def unnormalizedLogPdf(x: A): Double =
    dist.unnormalizedLogPdf(x)
  def gradLogPdf(x: A): Double
  def draw(): A = dist.draw
}

object GradDist {
  def gamma(dist: Gamma) = new GradDist[Double](dist) {
    def gradLogPdf(x: Double): Double =
      (dist.shape - 1.0) / x - 1.0 / dist.scale
  }

  def normal(dist: Gaussian) = new GradDist[Double](dist) {
    def gradLogPdf(x: Double): Double =
      - (x - dist.mu) / dist.sigma
  }
}
