package examples

import gp.core._
import breeze.stats.distributions._
import breeze.stats.mean
import breeze.linalg.{DenseVector, DenseMatrix, diag}
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
import math.{min, max}
import dlm.core.model._

trait TemperatureDlm {
  implicit val jodaDateTime: CellCodec[DateTime] = {
    val format = DateTimeFormat.forPattern("yyyy-MM-dd'T'HH:mm:ss'Z'")
    CellCodec.from(s => DecodeResult(format.parseDateTime(s)))(d => format.print(d))
  }

  case class Temperature(
    name:     String,
    date:     DateTime,
    obs:      Option[Double],
    lon:      Double,
    lat:      Double)

  // specify GP Parameters
  val initParams = GaussianProcess.Parameters(
    MeanParameters.zero,
    Vector(KernelParameters.se(1.0, 0.001), KernelParameters.white(1.0))
  )
  val covFn = KernelFunction.apply(initParams.kernelParameters)
  val dist = Location.euclidean _

  val test_sensor = "new_new_emote_2604"

  // read in raw data
  val rawData = Paths.get("data/daily_average_temp.csv")
  val reader = rawData.asCsvReader[Temperature](rfc.withHeader)
  val data: Vector[Temperature] = reader.
    collect {
      case Right(a) => a
    }.
    toVector

  val trainingData = data.
    filter(_.name != test_sensor)

  val testData = data.
    filter(_.name == test_sensor)

  // data formatted for DLM
  val ys = trainingData.
    groupBy(_.date).
    map { case (t, temps) =>
      Data(t.getMillis() / (60.0 * 60.0 * 1000.0 * 24), // time in days
           DenseVector(temps.map(s => s.obs).toArray))
    }.
    toVector.
    sortBy(_.time).
    zipWithIndex.
    map { case (d, i) => d.copy(time = i.toDouble) }


  implicit val randMonad = new Monad[Rand] {
    def pure[A](x: A): Rand[A] = Rand.always(x)
    def flatMap[A, B](fa: Rand[A])(f: A => Rand[B]): Rand[B] =
      fa flatMap f

    def tailRecM[A, B](a: A)(f: A => Rand[Either[A, B]]): Rand[B] =
      f(a).draw match {
        case Left(a1) => tailRecM(a1)(f)
        case Right(b) => Rand.always(b)
      }
  }

  val seasonalDlm = Dlm.polynomial(1) |+| Dlm.seasonal(7, 5)
  val mvDlmShareState = seasonalDlm.
    copy(f = (t: Double) => List.fill(8)(seasonalDlm.f(t)).
           reduce((a, b) => DenseMatrix.horzcat(a, b)))
  val model = DlmGp.Model(None, mvDlmShareState, dist)

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

  val file = s"data/temperature_dlm_share_state_0.csv"

  val dlmPs = Streaming.readCsv(file).
    map(_.map(_.toDouble).toList).
    drop(1000).
    map(x => DlmParameters.fromList(8, 11)(x)).
    via(Streaming.meanParameters(8, 11))
}

object FitTemperatureDlm extends App with TemperatureDlm {
  implicit val system = ActorSystem("temperature_dlm")
  implicit val materializer = ActorMaterializer()

  val priorV = InverseGamma(3.0, 3.0)
  val priorW = InverseGamma(3.0, 3.0)

  val prior = for {
    v <- Applicative[Rand].replicateA(8, priorV)
    w <- Applicative[Rand].replicateA(11, priorW)
    m0 = Array.fill(11)(0.0)
    c0 = Array.fill(11)(10.0)
  } yield DlmParameters(
      v = diag(DenseVector(v.toArray)),
      w = diag(DenseVector(w.toArray)),
      m0 = DenseVector(m0),
      c0 = diag(DenseVector(c0)))

  val iters = GibbsSampling.sampleSvd(mvDlmShareState, priorV,
                                   priorW, prior.draw, ys)

  def formatParameters(s: GibbsSampling.State) =
    s.p.toList

  Streaming.writeParallelChain(
    iters, 2, 100000, "data/temperature_dlm_share_state", formatParameters).
    runWith(Sink.onComplete { s =>
              println(s)
              system.terminate()
            })
}

object ForecastTemperatureDlm extends App with TemperatureDlm {
  implicit val system = ActorSystem("forecast_temperature")
  implicit val materializer = ActorMaterializer()

  val times: Vector[DateTime] = data.map(_.date)

  val out = new java.io.File(s"data/forecast_temperature.csv")
  val headers = rfc.withHeader(false)

  (for {
     p <- dlmPs.runWith(Sink.head)
     filtered = KalmanFilter.filterDlm(mvDlmShareState, ys, p)
     summary = filtered.flatMap(kf => (for {
                                         f <- kf.ft
                                         q <- kf.qt
                                       } yield Dlm.summariseForecast(0.75)(f, q)).
                                  map(f => f.flatten.toList))
     toWrite: Vector[(DateTime, List[Double])] = times.
     zip(summary).
     map(d => (d._1, d._2))
     io = out.writeCsv(toWrite, headers)
   } yield io).
    onComplete(_ => system.terminate())
}

object StateTemperatureDlm extends App with TemperatureDlm {
  implicit val system = ActorSystem("forecast_temperature")
  implicit val materializer = ActorMaterializer()

  val out = new java.io.File(s"data/state_temperature_dlm.csv")
  val sensorNames = trainingData.map(_.name).distinct.toList
  val headers = rfc.withHeader("day" :: sensorNames: _*)

  val res = for {
    p <- dlmPs.runWith(Sink.head)
    state = Smoothing.ffbsDlm(mvDlmShareState, ys, p).sample(1000).
      map(_.map(x => x.copy(sample = mvDlmShareState.f(x.time).t * x.sample)))
    times: Vector[DateTime] = trainingData.map(_.date)
    meanState = state.transpose.zip(times).
      map { case (x, t) => (t, x.map(_.sample).reduce(_ + _).map(_ / 1000).data.toList) }.toList
    io = out.writeCsv(meanState, headers)
  } yield io

  res.onComplete{ s =>
    println(s)
    system.terminate()
  }
}

object TemperatureDlmGp extends App with TemperatureDlm {
  implicit val system = ActorSystem("dlmgp_temperature_residuals")
  implicit val materializer = ActorMaterializer()

  def proposal(delta: Double)(ps: Vector[KernelParameters]) = ps.traverse(p =>
    p match {
      case SquaredExp(h, s) =>
        for {
          z1 <- Gaussian(0.0, delta)
          newh = h * math.exp(z1)
          z2 <- Gaussian(0.0, delta)
          newS = s * math.exp(z2)
        } yield KernelParameters.se(newh, newS)
      case White(s) =>
        for {
          z <- Gaussian(0.0, delta)
          news = s * math.exp(z)
        } yield KernelParameters.white(news)
    })

  def priorKernel(ps: Vector[KernelParameters]) = ps.map(p => p match {
    case SquaredExp(h, sigma) =>
      InverseGamma(3, 3).logPdf(h) +
      InverseGamma(3, 3).logPdf(sigma)
    case White(s) =>
      InverseGamma(3.0, 3).logPdf(s)
  }).sum

  // create a vector of time series grouped by location
  val trainingDataDlmGp: Vector[(Location[Double], Vector[Data])] = trainingData.
    groupBy(a => (a.lon, a.lat)).
    map { case (l, ts) =>
      (Two(l._1, l._2), ts.groupBy(_.date).
         map { case (t, d) => Data(t.getMillis / (60.0 * 60.0 * 1000.0 * 24.0),
                                   DenseVector(d.map(_.obs).toArray)) }.toVector)
    }.
    toVector

  val priorV = InverseGamma(3.0, 3.0)
  val priorW = InverseGamma(3.0, 3.0)

  val prior = for {
    v <- Applicative[Rand].replicateA(8, priorV) // a different measurement variance for each sensor
    w <- Applicative[Rand].replicateA(11, priorW)
    m0 <- Applicative[Rand].replicateA(11, Gaussian(0.0, 1.0))
    c0 <- Applicative[Rand].replicateA(11, InverseGamma(3.0, 3.0))
    h <- InverseGamma(3.0, 4.0)
    sigma <- InverseGamma(3.0, 3.0)
    sigmaY <- InverseGamma(3.0, 3.0)
  } yield DlmGp.Parameters(
    DlmParameters(
      v = diag(DenseVector(v.toArray)), // ignored
      w = diag(DenseVector(w.toArray)),
      m0 = DenseVector(m0.toArray),
      c0 = diag(DenseVector(c0.toArray))),
    GaussianProcess.Parameters(
      MeanParameters.zero,
      Vector(KernelParameters.se(h, sigma), KernelParameters.white(sigmaY))
    )
  )

  val iters = FitDlmGp.sample(priorV, priorW, priorKernel, proposal(0.05),
                              model, trainingDataDlmGp, prior.draw)

  def formatParameters(s: FitDlmGp.State): List[Double] = ???

  Streaming.writeParallelChain(
    iters, 2, 100000, "data/temperature_dlmgp_share_state", formatParameters).
    runWith(Sink.onComplete { s =>
              println(s)
              system.terminate()
            })
}

// Fit a zero mean covariance function to the residuals of the Temperature DLM
object FitTemperatureResiduals extends App {
  implicit val system = ActorSystem("gp_temperature_residuals")
  implicit val materializer = ActorMaterializer()

  implicit val randMonad = new Monad[Rand] {
    def pure[A](x: A): Rand[A] = Rand.always(x)
    def flatMap[A, B](fa: Rand[A])(f: A => Rand[B]): Rand[B] =
      fa flatMap f

    def tailRecM[A, B](a: A)(f: A => Rand[Either[A, B]]): Rand[B] =
      f(a).draw match {
        case Left(a1) => tailRecM(a1)(f)
        case Right(b) => Rand.always(b)
      }
  }

  implicit val jodaDateTime: CellCodec[DateTime] = {
    val format = DateTimeFormat.forPattern("yyyy-MM-dd'T'HH:mm:ss'Z'")
    CellCodec.from(s => DecodeResult(format.parseDateTime(s)))(d => format.print(d))
  }

  case class Residuals(
    day: DateTime,
    name: String,
    fitted: Double,
    temp: Double,
    lon: Double,
    lat: Double,
    residual: Double)

  val initParams = GaussianProcess.Parameters(
    MeanParameters.zero,
    Vector(KernelParameters.se(0.5, 1.0), KernelParameters.white(0.01))
  )
  val covFn = KernelFunction.apply(initParams.kernelParameters)
  val dist = Location.euclidean _

  // read in residual data and convert to a data point
  val rawData1 = Paths.get("data/dlm_temperature_residuals.csv")
  val reader1 = rawData1.asCsvReader[Residuals](rfc.withHeader)
  val residuals: Vector[GaussianProcess.Data] = reader1.
    collect {
      case Right(a) => a
    }.
    toVector.
    map(d => GaussianProcess.Data(Two(d.lon, d.lat), d.residual)).
    groupBy(_.x).
    map { case (l, ds) => GaussianProcess.Data(l, mean(ds.map(_.y))) }.
    toVector

  def proposal(delta: Double)(ps: Vector[KernelParameters]) = ps.traverse(p => p match {
    case SquaredExp(h, s) =>
      for {
        z1 <- Gaussian(0.0, delta)
        newh = h * math.exp(z1)
        z2 <- Gaussian(0.0, delta)
        newS = s * math.exp(z2)
      } yield KernelParameters.se(newh, newS)
    case White(s) =>
      for {
        z <- Gaussian(0.0, delta)
        news = s * math.exp(z)
      } yield KernelParameters.white(news)
  })

  // TODO: Think about more informative prior distributions
  val priorh = InverseGamma(10, 2)
  val priorSigma = InverseGamma(3, 5)
  val priorSigmaY = InverseGamma(10.0, 2.0)

  def priorKernel(ps: Vector[KernelParameters]) = ps.map(p => p match {
    case SquaredExp(h, sigma) =>
      priorh.logPdf(h) +
      priorSigma.logPdf(sigma)
    case White(s) =>
      priorSigmaY.logPdf(s)
  }).sum

  val prior = for {
    v <- priorSigmaY
    h <- priorh
    sigma <- priorSigma
  } yield GaussianProcess.Parameters(
    MeanParameters.zero,
    Vector(KernelParameters.se(h, sigma), KernelParameters.white(v))
  )

  // get iterations from command line argument
  val nIters: Int = args.lift(0).map(_.toInt).getOrElse(100000)

  val iters = Mcmc.sample(residuals, priorKernel,
    Gaussian(1.0, 1.0), proposal(0.05), dist, prior.draw)

  def formatParameters(p: GaussianProcess.Parameters): List[Double] = {
    p.kernelParameters.flatMap(_.toList).toList ::: p.meanParameters.toList
  }

  Streaming.writeParallelChain(
    iters, 2, nIters, "data/temperature_gp_residuals", formatParameters).
    runWith(Sink.onComplete { s =>
              println(s)
              system.terminate()
            })
}

// perform a forecast using the DLM and GP
// currently each training sensor has it's own measurement variance
// hence predicting measurements at a new sensor location is impossible without
// specifying the measurement variance -- this could be learned online using a conjugate filter...
object ForecastTestDlm extends App with TemperatureDlm {
  implicit val system = ActorSystem("forecast_dlm_gp")
  implicit val materializer = ActorMaterializer()

  // extract the times from the test data
  val times: Vector[DateTime] = testData.map(_.date)

  // read test location
  val locationFile = Paths.get("data/test_temperature_location.csv")
  val locReader = locationFile.asCsvReader[(String, Double, Double)](rfc.withHeader)
  val testLocation = locReader.
    collect {
      case Right(a) => Two(a._2, a._3)
    }.
    toVector

  val trainingDataDlmGp: Vector[Data] = trainingData.
    groupBy(_.date).
    map { case (date, ts) =>
      Data(date.getMillis / (24 * 60 * 60 * 1000), DenseVector(ts.map(_.obs).toArray)) }.
    toVector

  val trainingLocations = trainingData.
    map(t => Two(t.lon, t.lat)).
    distinct

  val testDataDlmGp = testData.
    map(d => Data(d.date.getMillis / (24 * 60 * 60 * 1000), DenseVector(d.obs)))

  val out = new java.io.File("data/temperature_test_predictions_dlmgp.csv")
  val headers = rfc.withHeader("day", "mean", "lower", "upper")

  (for {
     // get mean parameters
     dlmP <- dlmPs.runWith(Sink.head)
     p = DlmGp.Parameters(dlmP, gpParams)

     // perform prediction at test location for each day
     fitted = DlmGp.forecast(trainingDataDlmGp, trainingLocations, testDataDlmGp,
                             testLocation, model, seasonalDlm, p).
     map(a => for {
           m <- a.mean
           c <- a.cov
         } yield Gaussian(m.data.head, c.data.head)).
     flatten

     _ = fitted.foreach(println)

    predictions = Predict.predict(fitted, 0.75)

     // write out predictions
     io = out.writeCsv(times zip predictions, headers)
   } yield io).
    onComplete{
      s => println(s)
      system.terminate()
    }
}

// perform kriging on a regular grid to interpolate the temperature
// could use more temperature sensors for this too
object KrigTemperature extends App with TemperatureDlm {
  implicit val system = ActorSystem("forecast_dlm_gp")
  implicit val materializer = ActorMaterializer()

  val times: Vector[DateTime] = trainingData.map(_.date)

  // 1. Get a grid of locations
  // 2. Use DlmGp.forecast to forecast at every location
  // 3. Write the mean and covariance for each

  val trainingLocations = trainingData.
    map(t => Two(t.lon, t.lat)).
    distinct

  val lons: Vector[Double] = trainingLocations.map(_.x)
  val lats: Vector[Double] = trainingLocations.map(_.y)

  val locations = DlmGp.getGridLocations(lons, lats, 0.005)

  locations foreach println

  val trainingDataDlmGp: Vector[Data] = trainingData.
    groupBy(_.date).
    map { case (date, ts) =>
      Data(date.getMillis / (24 * 60 * 60 * 1000), DenseVector(ts.map(_.obs).toArray)) }.
    toVector

  val out = new java.io.File("data/temperature_kriging_dlmgp.csv")
  val headers = rfc.withHeader("day", "lon", "lat", "mean")

  (for {
     // get mean parameters
     dlmP <- dlmPs.runWith(Sink.head)
     p = DlmGp.Parameters(dlmP, gpParams)

     // perform prediction at test location for each day
     fitted = DlmGp.forecastLocations(trainingDataDlmGp, trainingLocations,
                             locations, model, seasonalDlm, p).
     map(_._2)

     // predictions = Predict.predict(fitted, 0.75)

     // write out predictions
     locsPred: Vector[(DateTime, Double, Double, Double)] = (locations zip fitted).
       zipWithIndex.
       flatMap { case ((Two(x, y), f), i) => f.data.toVector.map(fi => (times(i), x, y, fi)) }

     io = out.writeCsv(locsPred, headers)
   } yield io).
    onComplete{
      s => println(s)
      system.terminate()
    }
}
