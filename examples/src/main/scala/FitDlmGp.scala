package examples

import breeze.stats.distributions._
import breeze.linalg.{DenseVector, DenseMatrix, diag}
import core.dlm.model._
import gp.core._
import cats.implicits._

/**
  * Fit a DLM GP to data
  * The parameters of the DLM (W, m0, c0 and the state x_{0:T}) should
  * be shared amongst all of the sensors. Then the spatial differences
  * will be accounted for by the Gaussian Process
  */
object FitDlmGp {
  case class State(
      p: DlmGp.Parameters,
      state: Vector[(Double, DenseVector[Double])]
  )

  /**
    * Calculate Y_t - F(X_t)
    * @param ys a vector of observations
    * @param state a vector of the DLM state X_t
    */
  def residual(
      ys: Vector[DlmGp.Data],
      state: Vector[(Double, DenseVector[Double])]
  ): Vector[(Double, Vector[GaussianProcess.Data])] = {

    ys.map(d =>
      (d.time, for {
        (loc, measurement) <- d.locations zip d.measurement.data
      } yield GaussianProcess.Data(loc, measurement)))
  }

  /**
    * Convert a DLM GP to a Time series of mean values
    * @param ys a vector of observations from m locations
    * @return a time series vector containing the mean value across all
    * locations
    */
  def gptoTimeSeries(ys: Vector[DlmGp.Data]) = {

    ys.map { d =>
      Dlm.Data(d.time, d.measurement.map(_.some))
    }
  }

  /**
    *
    */
  def dlmParameters(kxx: DenseMatrix[Double],
                    p: DlmGp.Parameters): DlmParameters = {

    DlmParameters(kxx, p.w, p.m0, p.c0)
  }

  def step(priorW: InverseGamma,
           priorKern: Vector[KernelParameters] => Double,
           propKern: Vector[KernelParameters] => Rand[Vector[KernelParameters]],
           mod: DlmGp.Model,
           ys: Vector[DlmGp.Data],
           xs: Vector[Location[Double]]) = { st: State =>
    val covFn = KernelFunction(st.p.gp.kernelParameters)
    val nugget = diag(DenseVector.fill(xs.size)(1e-3))
    val kxx = KernelFunction.buildCov(xs, covFn, mod.dist) + nugget
    val dlmP = dlmParameters(kxx, st.p)

    val dlm = toDlmModel(mod)

    for {
      state <- Smoothing.ffbsDlm(dlm, gptoTimeSeries(ys), dlmP)
      system <- GibbsSampling.sampleSystemMatrix(priorW, state, dlm.g)
      res = residual(ys, state)
      kp <- KernelParameters.sample(res.head._2, mod.dist, priorKern, propKern)(
        st.p.gp)
    } yield
      State(
        st.p.copy(w = system,
                  gp = st.p.gp.copy(kernelParameters = kp.kernelParameters)),
        state)
  }

  /**
    * The DLMs share the same state, change the observation matrix to reflect this
    */
  def toDlmModel(mod: DlmGp.Model): DlmModel = {
    val newf = mod.dim
      .map { p => (time: Double) =>
        List
          .fill(p)(mod.dlm.f(time))
          .reduce((a, b) => DenseMatrix.horzcat(a, b))
      }
      .getOrElse((time: Double) => mod.dlm.f(time))

    mod.dlm.copy(f = newf)
  }

  def sample(
      priorW: InverseGamma,
      priorKern: Vector[KernelParameters] => Double,
      propKern: Vector[KernelParameters] => Rand[Vector[KernelParameters]],
      mod: DlmGp.Model,
      ys: Vector[DlmGp.Data],
      xs: Vector[Location[Double]],
      initParams: DlmGp.Parameters) = {

    val covFn = KernelFunction.apply(initParams.gp.kernelParameters)
    val nugget = diag(DenseVector.fill(xs.size)(1e-3))
    val kxx = KernelFunction.buildCov(xs, covFn, mod.dist) + nugget

    val dlmP = dlmParameters(kxx, initParams)
    val dlmMod = toDlmModel(mod)

    val init = for {
      st <- Smoothing.ffbsDlm(dlmMod, gptoTimeSeries(ys), dlmP)
    } yield State(initParams, st)

    MarkovChain(init.draw)(step(priorW, priorKern, propKern, mod, ys, xs))
  }
}
