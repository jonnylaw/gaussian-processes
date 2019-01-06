package com.github.jonnylaw.gp.examples

import breeze.stats.distributions._
import breeze.linalg.{DenseVector, DenseMatrix, diag}
import dlm.core.model._
import com.github.jonnylaw.gp._
import cats.implicits._

/**
  * Fit a DLM GP to data
  * The parameters of the DLM (W, m0, c0 and the state x_{0:T}) should
  * be shared amongst all of the sensors. Then the spatial differences
  * will be accounted for by the Gaussian Process
  */
object FitDlmGp {
  case class State(
    p:   DlmGp.Parameters,
    dlm: GibbsSampling.State)

  /**
    * Calculate Y_t - F.t * X_t and return a collection of locations and
    * measurements
    * @param ys a collection of time series measured at different locations
    * vector of observations
    * @param state a vector of the DLM state X_t
    * @param f the time dependent observation matrix
    * @return a vector of
    */
  def residual(
    ys: Vector[(Location[Double], Vector[Data])],
    state: Vector[SamplingState],
    f: Double => DenseMatrix[Double]): Vector[GaussianProcess.Data] = {

    val locations = ys.map(_._1)
    val ts: Vector[Data] = combineTs(ys)
    val residual: Vector[Array[Option[Double]]] = (ts.map(_.observation) zip state.tail).
      map { case (d, s) =>
        val ft: DenseVector[Double] = f(s.time).t * s.sample
        d.data.zipWithIndex.map {
          case (Some(y), i) => Some(y - ft(i))
          case _ => None
        }
      }

    residual.flatMap(r => r.zipWithIndex.
      map { case (Some(r), i) => GaussianProcess.Data(locations(i), r) })
  }

  /**
    * Combine time series at multiple locations into one multivariate time
    * series
    */
  def combineTs(ys: Vector[(Location[Double], Vector[Data])]): Vector[Data] =
    ys.flatMap(_._2).groupBy(_.time).
      map { case (t, ds) => Data(t, DenseVector(ds.flatMap(_.observation.data).toArray)) }.toVector

  def step(
    priorV:    InverseGamma,
    priorW:    InverseGamma,
    priorKern: Vector[KernelParameters] => Double,
    propKern:  Vector[KernelParameters] => Rand[Vector[KernelParameters]],
    mod:       DlmGp.Model,
    ys:        Vector[(Location[Double], Vector[Data])]) = { st: State =>

    val covFn = KernelFunction(st.p.gp.kernelParameters)
    val xs = ys.map(_._1)
    val nugget = diag(DenseVector.fill(xs.size)(1e-3))
    val kxx = KernelFunction.buildCov(xs, covFn, mod.dist) + nugget
    val ts = combineTs(ys)

    for {
      dlmSt <- GibbsSampling.stepSvd(mod.dlm, priorV, priorW, ts)(st.dlm)
      res = residual(ys, dlmSt.state, mod.dlm.f)
      kp <- KernelParameters.sample(res, mod.dist, priorKern, propKern)(st.p.gp)
    } yield State(st.p.copy(dlm = dlmSt.p.copy(v = kxx),
                            gp = st.p.gp.copy(kernelParameters = kp.kernelParameters)),
      dlmSt)
  }

  /**
    * The DLMs share the same state,
    * change the observation matrix to reflect this
    */
  def toDlm(mod: DlmGp.Model): Dlm = {
    val newf = mod.dim.
      map { p =>
        (time: Double) => List.fill(p)(mod.dlm.f(time)).
          reduce((a, b) => DenseMatrix.horzcat(a, b))
      }.
      getOrElse((time: Double) => mod.dlm.f(time))
    mod.dlm.copy(f = newf)
  }

  /**
    * Perform metropolis-with-Gibbs sampling to determine the posterior
    * distribution of the parameters of a DLM GP
    * @param priorV
    * @param priorW
    * @param priorKern
    * @param propKern the proposal distribution for the kernel parameters
    * @param mod a DLM GP model
    * @param ys
    */
  def sample(
    priorV:     InverseGamma,
    priorW:     InverseGamma,
    priorKern:  Vector[KernelParameters] => Double,
    propKern:   Vector[KernelParameters] => Rand[Vector[KernelParameters]],
    mod:        DlmGp.Model,
    ys:         Vector[(Location[Double], Vector[Data])],
    initParams: DlmGp.Parameters) = {

    val covFn = KernelFunction.apply(initParams.gp.kernelParameters)
    val xs = ys.map(_._1)
    val nugget = diag(DenseVector.fill(xs.size)(1e-3))
    val kxx = KernelFunction.buildCov(xs, covFn, mod.dist) + nugget

    val dlmMod = toDlm(mod)
    val ts = combineTs(ys)

    val init = for {
      st <- SvdSampler.ffbsDlm(dlmMod, ts, initParams.dlm.copy(v = kxx))
      _ = st.foreach(println)
      initSt = GibbsSampling.State(initParams.dlm.copy(v = kxx), st)
    } yield State(initParams, initSt)

    MarkovChain(init.draw)(step(priorV, priorW, priorKern,
                                propKern, mod, ys))
  }
}

