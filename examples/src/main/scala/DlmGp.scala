package examples

import core.dlm.model._
import gp.core._
import breeze.linalg.{DenseMatrix, DenseVector}
import breeze.stats.distributions._

object DlmGp {

  /**
    * A model from the DLM GP
    * @param dim the dimension of the DLM GP, if this is omitted a composed DLM
    * will be used
    * @param dlm a DLM model (either composed or a single model)
    * @param dist a distance function for the Gaussian Process Model
    */
  case class Model(dim: Option[Int],
                   dlm: DlmModel,
                   dist: (Location[Double], Location[Double]) => Double)

  /**
    * Parameters for the DLM GP
    * @param v the covaraiance of the system noise matrix
    * @param gp parameters of the gaussian process
    */
  case class Parameters(w: DenseMatrix[Double],
                        m0: DenseVector[Double],
                        c0: DenseMatrix[Double],
                        gp: GaussianProcess.Parameters)

  case class Data(time: Double,
                  locations: Vector[Location[Double]],
                  measurement: DenseVector[Double])

  /**
    * Simulate a single step of the DLM GP
    * @param mod a DLM GP model where the observation covariance is modelled with
    * a gaussian process
    * @param params the parameters of the DLM GP
    * @param xs the unique locations to simulate at
    * @param time the current time of the observation
    * @param theta a realisation of the state vector at time t-1
    */
  def simStep(mod: Model, params: Parameters, xs: Vector[Location[Double]])(
      time: Double,
      theta: DenseVector[Double]) = {

    val f = mod.dim
      .map { p =>
        List
          .fill(p)(mod.dlm.f(time))
          .reduce((a, b) => DenseMatrix.horzcat(a, b))
      }
      .getOrElse(mod.dlm.f(time))

    for {
      v <- GaussianProcess.draw(xs, mod.dist, params.gp)
      w <- MultivariateGaussianSvd(DenseVector.zeros[Double](params.w.cols),
                                   params.w)
      theta1 = mod.dlm.g(1.0) * theta + w
      y = f.t * theta1 + v
    } yield (Data(time, xs, y), theta1)
  }

  /**
    *
    */
  def simulate(mod: Model, p: Parameters, xs: Vector[Location[Double]]) = {

    val x0 = MultivariateGaussian(p.m0, p.c0).draw
    val init = (Data(0.0, xs, DenseVector.zeros[Double](xs.size)), x0)

    MarkovChain(init) {
      case (d, theta) => simStep(mod, p, xs)(d.time + 1.0, theta)
    }
  }
}
