package com.github.jonnylaw.gp

import cats._

sealed trait Location[+A] {
  def toVector: Vector[A]
}
case class One[A](x: A) extends Location[A] {
  def toVector: Vector[A] = Vector(x)
}
case class Two[A](x: A, y: A) extends Location[A] {
  def toVector: Vector[A] = Vector(x, y)
}

object Location {
  implicit val tol = 1e-3

  implicit val eq = new Eq[Location[Double]] {
    def eqv(a: Location[Double], b: Location[Double]) = (a, b) match {
      case (One(x1), One(x2)) => math.abs(x1 - x2) < tol
      case (Two(x1, y1), Two(x2, y2)) =>
        math.abs(x1 - x2) < tol && math.abs(y1 - y2) < tol
      case _ => false
    }
  }

  def euclidean(a: Location[Double], b: Location[Double]): Double =
    (a, b) match {
      case (One(x1), One(x2)) => math.sqrt(math.pow(x1 - x2, 2))
      case (Two(x1, y1), Two(x2, y2)) =>
        math.sqrt(math.pow(x1 - x2, 2) + math.pow(y1 - y2, 2))
      case _ => throw new Exception("Dimension mismatch")
    }
}
