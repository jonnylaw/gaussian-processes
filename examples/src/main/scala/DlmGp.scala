package com.github.jonnylaw.gp.examples

import dlm.core.model._
import com.github.jonnylaw.gp._
import breeze.linalg.{DenseMatrix, DenseVector, diag, svd, inv}
import breeze.stats.distributions._
import cats.implicits._

object DlmGp {

  /**
    * A model from the DLM GP
    * @param dim the dimension of the DLM GP, if this is omitted a composed DLM
    * will be used?
    * @param dlm a DLM model (either composed or a single model)
    * @param dist a distance function for the Gaussian Process Model
    */
  case class Model(
    dim:  Option[Int],
    dlm:  Dlm,
    dist: (Location[Double], Location[Double]) => Double)

  /**
    * Parameters for the DLM GP
    * @param v the covaraiance of the system noise matrix
    * @param gp parameters of the gaussian process
    */
  case class Parameters(
    dlm: DlmParameters,
    gp: GaussianProcess.Parameters)

  case class DlmGpData(
    time:         Double,
    locations:    Vector[Location[Double]],
    measurement:  DenseVector[Double])

  /**
    * Simulate a single step of the DLM GP
    * @param mod a DLM GP model where the observation covariance is modelled with
    * a gaussian process
    * @param params the parameters of the DLM GP
    * @param xs the unique locations to simulate at
    * @param time the current time of the observation
    * @param theta a realisation of the state vector at time t-1
    */
  def simStep(
    mod:    Model,
    params: Parameters,
    xs:     Vector[Location[Double]])
    (time: Double, theta: DenseVector[Double]) = ???

  //   val f = mod.dim.
  //     map { p =>
  //       List.fill(p)(mod.dlm.f(time)).reduce((a, b) => DenseMatrix.horzcat(a, b))
  //     }.
  //     getOrElse(mod.dlm.f(time))

  //   for {
  //     ys <- Dlm.simStep(mod.dlm, params.dlm)(theta, time, 1.0)
  //     vs <- GaussianProcess.draw(xs, mod.dist, params.gp)
  //     // TODO: Add GP to DLM
  //   } yield (DlmGpData(ys._1.time, xs, ys._1.observation.map(_ + vs)), ys._2)
  // }

  /**
    *
    */
  def simulate(mod: Model, p: Parameters, xs: Vector[Location[Double]]) = {

    val x0 = MultivariateGaussian(p.dlm.m0, p.dlm.c0).draw
    val init = (DlmGpData(0.0, xs, DenseVector.zeros[Double](xs.size)), x0)

    MarkovChain(init) {
      case (d, theta) => simStep(mod, p, xs)(d.time + 1.0, theta)
    }
  }

  case class State(
    kfTraining: KfState,
    kfTest: InverseGammaState,
    mean: Option[DenseVector[Double]],
    cov: Option[DenseMatrix[Double]]
  )

  /**
    * Take away the one-step predictive mean from the measurements
    * @param data a time series object with encoded missing data
    * @param kf the state of a Kalman Filter object
    * @return the residuals
    */
  def residuals(
    data: Data,
    ft: DenseVector[Double]): Data = {

    data.copy(observation =
                DenseVector((data.observation.data zip ft.data).
                  map { case (yo, f) => yo.map(y => y - f) } ))
  }

  def forecastStep(
    testLocation: Vector[Location[Double]],
    trainingLocations: Vector[Location[Double]],
    mod: Model,
    testMod: Dlm,
    p: Parameters)(state: State, observed: (Data, Data)) = {

    val trainingDataTs = observed._2

    // perform one-step forecast for training data
    val kfTraining = KalmanFilter(KalmanFilter.advanceState(p.dlm, mod.dlm.g)).
      step(mod.dlm, p.dlm)(state.kfTraining, trainingDataTs)

    // take one-step predictive mean away from the training data
    val trainingResiduals = residuals(trainingDataTs, kfTraining.ft.get)

    // predict GP at test location given the training data
    val trainingData = toGpData(trainingLocations, trainingResiduals)
    val gp = Predict.fit(testLocation, trainingData, mod.dist, p.gp)

    // perform one-step forecast for test data using the conjugate filter
    val testData = observed._1
    val priorV = InverseGamma(3.0, 2.0)
    val kfTest = ConjugateFilter(priorV, ConjugateFilter.advanceState(p.dlm, testMod.g)).
      step(testMod, p.dlm)(state.kfTest, testData)

    val mean = kfTest.kfState.ft.map (f => f + DenseVector(gp.map(_._2.mean).toArray))
    val cov = kfTest.kfState.qt.map(q => q + diag(DenseVector(gp.map(_._2.variance).toArray)))

    // println(s"mean $mean")
    // println(s"cov $cov")

    State(kfTraining, kfTest, mean, cov)
  }

  /**
    * Convert a time series with encoded missing data into a
    * times series with location information
    */
  def toGpData(
    locations: Vector[Location[Double]],
    data: Data): Vector[GaussianProcess.Data] = {

    val nonMissing = KalmanFilter.indexNonMissing(data.observation)
    val ls = nonMissing map (locations(_))

    for {
      (y, l) <- data.observation.data.flatten.toVector zip ls
    } yield GaussianProcess.Data(l, y)
  }

  /**
    * Perform a forecast at a new location for the DLM GP
    * by calculating the mean of the one-step forecast for the dlm f_t = F_ta_t
    * then adding the mean of the GP prediction at that location for that given
    * time point
    * @param trainingData the locations and associated time series
    * for the Gaussian Process training data
    * @param testData the test time series data
    * @param testLocation
    */
  def forecast(
    trainingData: Vector[Data],
    trainingLocations: Vector[Location[Double]],
    testData: Vector[Data],
    testLocation: Vector[Location[Double]],
    mod: Model,
    testMod: Dlm,
    p: Parameters) = {

    val priorV = InverseGamma(3.0, 4.0)
    val t0 = testData.map(_.time).reduceLeftOption((t0, d) => math.min(t0, d))

    // val root = svd(p.dlm.c0)
    // val dc0 = root.singularValues.map(math.sqrt)
    // val uc0 = root.rightVectors.t

    val ft = mod.dlm.f(t0.get).t * p.dlm.m0
    // val kf = SvdState(t0.get - 1, p.dlm.m0, dc0, uc0, p.dlm.m0, dc0, uc0, ft)
    val kf = KfState(t0.get - 1, p.dlm.m0, p.dlm.c0, p.dlm.m0, p.dlm.c0, None, None, 0.0)
    val v0 = Vector.fill(testLocation.size)(priorV)
    val ckf = InverseGammaState(KfState(t0.get - 1.0, p.dlm.m0, p.dlm.c0,
                                        p.dlm.m0, p.dlm.c0, None, None, 0.0), v0)
    val init = State(kf, ckf, None, None)

    (testData zip trainingData).
      scanLeft(init)(forecastStep(testLocation, trainingLocations,
                                  mod, testMod, p))
  }

  /**
    * Calculate the minimum and maximum of a sequence with a single fold
    */
  def minMax(m: Seq[Double]): (Double, Double) =
    m.tail.foldLeft((m(0), m(0))){ case ((mi, ma), a) => (math.min(mi, a), math.max(ma, a)) }

  /**
    * Create a grid of locations within a bounding box
    * @param lons a collection of longitudes
    * @param lats a collection of latitudes
    * @param increment the step on the grid of latitude and logitude
    * @return
    */
  def getGridLocations(
    lons: Vector[Double],
    lats: Vector[Double],
    increment: Double): Vector[Location[Double]] = {

    val (minLon: Double, maxLon: Double) = minMax(lons)
    val (minLat: Double, maxLat: Double) = minMax(lats)

    (for {
       lo <- minLon until maxLon by increment
       la <- minLat until maxLat by increment
    } yield Two(lo, la)).toVector
  }

  /*
   * Calculate the one-step predictive variance for a
   * Kalman Filter implemented
   */
  def qtSvd(
    mod: Dlm,
    st: SvdState,
    sqrtVInv: DenseMatrix[Double]) = {
    val ft = mod.f(st.time)
    val rt = st.ur * diag(st.dr) * st.ur.t
    val v = inv(sqrtVInv.t * sqrtVInv)
    ft.t * rt * ft + v
  }
}
