package com.github.jonnylaw.gp

import breeze.linalg._
import breeze.numerics._

object Classify {
  /**
    * Multi class classification function returning the probability of fx
    * being in class C
    * @param fx an input evaluation of length c
    * @param c the class to consider
    * @return a Double representing the p(y = c | fx)
    */
  def softmax(fx: DenseVector[Double], c: Int): Double =
    exp(fx(c)) / sum(exp(fx))

  /**
    * Determine if an observation y belongs to class c
    */
  def delta(y: DenseVector[Int], c: Int): Int =
    y(c - 1)

  /**
    * Log-likelihood for multi-class classification
    * @param fxs a sequence of fitted functions
    * @param ys a sequence of one-hot encoded vectors
    * @param classes the total number of classes
    */
  def softmaxLl(
    fxs: Seq[Double],
    ys: Seq[DenseVector[Int]],
    classes: Int): Double = {

    val ls = for {
      (x, y) <- fxs zip ys
      c <- 0 until classes
      ll = delta(y, c) * (x - log(sum(fxs)))
    } yield ll

    ls.sum
  }

  /**
    * Encode the labels as "one-hot vectors" of length classes
    */
  def encodeLabels(y: Seq[Int], classes: Int): Seq[DenseVector[Double]] = {
    y.map { yi =>
      val m = DenseVector.zeros[Double](classes)
      m(yi) = 1.0
      m
    }
  }

  def objective(
    a: DenseVector[Double],
    fxs: Seq[DenseVector[Double]],
    ys: Seq[Int],
    classes: Int): Double = {
    val f = DenseVector.vertcat(fxs: _*)
    val y = DenseVector.vertcat(encodeLabels(ys, classes): _*)
      -0.5 * (a dot f) + (y dot f) + sum(
        log(fxs.map(fx => exp(fx)).reduce(_ + _)))
  }

  /**
    * Build a block diagonal matrix by combining two matrices of the same size
    */
  def blockDiagonal(
    a: DenseMatrix[Double],
    b: DenseMatrix[Double]): DenseMatrix[Double] = {

    val right = DenseMatrix.zeros[Double](a.rows, b.cols)
    val left = DenseMatrix.zeros[Double](b.rows, a.cols)

    DenseMatrix.vertcat(
      DenseMatrix.horzcat(a, right),
      DenseMatrix.horzcat(left, b)
    )
  }

  /**
    * Perform Newton-Raphson to determine the (unique) maximum value for f
    * GPML Algorithm 3.3
    * @param ys a sequence of labels of length n
    * @param ks a sequence of covariance matrices of length c
    * @param tol the tolerance of the Newton Method
    * @param fxs the values of the latent function for each class
    */
  def fit(
    ys: Seq[Int],
    ks: Seq[DenseMatrix[Double]],
    tol: Double,
    classes: Int)(fxs: Seq[DenseVector[Double]]) = {

    // build a list of classification matrices
    val pis: Seq[DenseVector[Double]] =
      fxs.map(fx => DenseVector(ys.map(y => softmax(fx, y)).toArray))

    // vertical concatenate the pis
    val pi: DenseMatrix[Double] = DenseMatrix.vertcat(pis.map(x => diag(x)): _*)
    val n = ys.size

    def loop(f: Seq[DenseVector[Double]],
             obj: Double): (Double, Seq[DenseVector[Double]]) = {
      val es = for {
        c <- 0 until classes
        d = diag(sqrt(pis(c)))
        k = ks(c)
        l = cholesky(DenseMatrix.eye[Double](n) + d * k * d)
        e = d * l.t \ (l \ d)
        z = sum(log(l))
      } yield (e, z)

      val d = diag(pis.reduce((a, b) => DenseVector.vertcat(a, b)))
      val m = cholesky(es.map(_._1).reduce(_ + _))
      val fs = DenseVector.vertcat(fxs: _*)
      val y = DenseVector.vertcat(encodeLabels(ys, classes): _*)
      val b = (d - pi * pi.t) * fs + y - DenseVector.vertcat(pis: _*)
      val e = es.map(_._1).reduce(blockDiagonal)
      val k = ks.reduce(blockDiagonal)
      val c = e * k * b
      val r =
        DenseMatrix.vertcat(Seq.fill(classes)(DenseMatrix.eye[Double](n)): _*)
      val a = b - c + e * r * m.t \ (m \ (r.t * c))
      val zs = es.map(_._2)
      val ll = -0.5 * (a dot fs) + (y dot fs) +
        sum(log(f.map(fx => exp(fx)).reduce(_ + _))) - sum(zs)

      val newF = ks.map(kc => kc * a)
      val newObj = objective(a, newF, ys, classes)

      if (abs(newObj - obj) < tol) {
        (ll, newF)
      } else {
        loop(newF, newObj)
      }
    }

    loop(fxs, 0.0)
  }

  /**
    * Perform prediction for a new unlabelled data point
    * @param k a sequence of independent covariance matrices length C
    * @param fxs the posterior mode from a training dataset calculated using the
    * Newton-Raphson Method using a Laplace Approximation
    * @param xs a new point to classify
    * @return a prediction vector of length C
    */
  def predict(k: Seq[DenseMatrix[Double]],
              fxs: Seq[DenseVector[Double]],
              xs: DenseVector[Double]): DenseVector[Double] = ???
}
