package com.github.jonnylaw.gp.examples

import com.github.jonnylaw.gp._
import breeze.stats.distributions._
import breeze.linalg.{DenseVector, diag}
import java.nio.file.Paths
import com.github.nscala_time.time.Imports._
import cats._
import cats.implicits._
import kantan.csv._
import kantan.csv.ops._
import kantan.csv.generic._
import akka.actor.ActorSystem
import akka.stream._
import scaladsl._
import scala.concurrent.ExecutionContext.Implicits.global
import dlm.core.model._

// perform kriging on a regular grid to interpolate the temperature
// could use more temperature sensors for this too
object KrigTemperature extends App with TemperatureDlm {
  val times: Vector[DateTime] = trainingData.map(_.date)

  // 1. Get a grid of locations
  // 2. Use DlmGp.forecast to forecast at every location
  // 3. Write the mean and covariance for each

  val trainingLocations = trainingData.
    map(t => Two(t.lon, t.lat)).
    distinct

  val lons: Vector[Double] = trainingLocations.map(_.x)
  val lats: Vector[Double] = trainingLocations.map(_.y)

  val locations = DlmGp.getGridLocations(lons, lats, 0.001)

  val paramsFile = Paths.get("data/temperature_gp_residuals_0.csv")
  val paramsReader = paramsFile.asCsvReader[List[Double]](rfc.withHeader)
  val chain = paramsReader.
    collect {
      case Right(a) => GaussianProcess.Parameters(
        MeanParameters.zero,
        Vector(KernelParameters.se(a(0), a(1)), KernelParameters.white(a(2)))
      )
    }.
    drop(1000).
    toStream.
    zipWithIndex.
    filter { case (_, i) => i % 20 == 0 }.
    map(_._1)

  // calculate mean of gpParameters parameters
  val gpParams: GaussianProcess.Parameters = chain.
    reduce { (a, b) =>
      GaussianProcess.Parameters(
        a.meanParameters add b.meanParameters,
        (a.kernelParameters zip b.kernelParameters).
          map { case (x, y) => x add y })
    }.
    map(p => p / chain.size)

  case class Residuals(
    date: DateTime,
    name: String,
    fitted: Double,
    temp: Double,
    lon: Double,
    lat: Double,
    residual: Double)

  val rawData1 = Paths.get("data/dlm_temperature_residuals.csv")
  val reader1 = rawData1.asCsvReader[Residuals](rfc.withHeader)
  val residuals: Vector[(DateTime, Vector[GaussianProcess.Data])] = reader1.
    collect {
      case Right(a) => a
    }.
    toVector.
    groupBy(_.date).
    map { case (date, ra) =>
      (date, ra.map(d => GaussianProcess.Data(Two(d.lon, d.lat), d.residual))) }.
    toVector.
    sortBy(_._1)

  def predict(p: GaussianProcess.Parameters): Vector[(Location[Double], Gaussian)] = {
    val meanResidual = residuals.
      flatMap(_._2).
      groupBy(_.x).
      map { case (location, ys) => GaussianProcess.Data(location, ys(30).y) }.
      toVector

    Predict.fit(locations, meanResidual, dist, gpParams)
  }

  // sample from the parameters and perform a fit
  val indices = Vector.fill(1000)(scala.util.Random.nextInt(chain.size))
  val parameters = indices.toVector.map(i => chain(i))
  parameters.take(10).foreach(println)
  val res: Vector[List[Double]] = parameters.
    flatMap(predict).
    groupBy(_._1).
    map { case (Two(lon, lat), ls) =>
      List(lon, lat, breeze.stats.mean(ls.map { case (l, g) => g.draw }), breeze.stats.variance(ls.map { case (l, g) => g.draw })) }.
    toVector

  val out = new java.io.File("data/temperature_kriging.csv")
  val headers = rfc.withHeader("lon", "lat", "mean", "variance")
  out.writeCsv(res, headers)
}
