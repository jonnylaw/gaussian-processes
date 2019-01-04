package com.github.jonnylaw.gp

import breeze.linalg.DenseVector

object MeanFunction {

  def plane(params: DenseVector[Double])(loc: Location[Double]) = loc match {
    case One(x)    => params(0) + x * params(1)
    case Two(x, y) => params(0) + x * params(1) + y * params(2)
  }

  def zero(loc: Location[Double]): Double = 0.0

  def apply(p: MeanParameters): Location[Double] => Double = p match {
    case Plane(beta) => plane(beta) _
    case Zero        => zero _
  }
}
