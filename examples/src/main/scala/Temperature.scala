package examples

import gp.core._
import breeze.stats.distributions._
import breeze.linalg.{DenseVector, diag}
import java.nio.file.Paths
import com.github.nscala_time.time.Imports._
import cats._
import cats.implicits._
import kantan.csv._
import kantan.csv.ops._
import kantan.csv.generic._
import dlm.core.model._

trait TemperatureModel {
  implicit val jodaDateTime: CellCodec[DateTime] = {
    val format = DateTimeFormat.forPattern("yyyy-MM-dd'T'HH:mm:ss'Z'")
    CellCodec.from(s => DecodeResult(format.parseDateTime(s)))(d =>
      format.print(d))
  }

  case class Temperature(
    name:     String,
    date:     DateTime,
    obs:      Option[Double],
    lon:      Double,
    lat:      Double)

  val initParams = GaussianProcess.Parameters(
    MeanParameters.plane(DenseVector(0.0, 0.0, 0.0)),
    Vector(KernelParameters.se(1.0, 0.001), KernelParameters.white(1.0))
  )
  val covFn = KernelFunction.apply(initParams.kernelParameters)
  val dist = Location.euclidean _

  val rawData = Paths.get("data/daily_average_temp.csv")
  val reader = rawData.asCsvReader[Temperature](rfc.withHeader)
  val data: Vector[Temperature] = reader.
    collect {
      case Right(a) => a
    }.
    toVector

  val test_sensor = "new_new_emote_2604"

  val testData = data.
    filter(_.name == test_sensor)

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
}

// TODO: Should I fit the mean function first
// or incorporate it into the Gibs sampler
object Temperature extends App with TemperatureModel {
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

  def priorKernel(ps: Vector[KernelParameters]) = ps.map(p => p match {
    case SquaredExp(h, sigma) =>
      InverseGamma(3, 5).logPdf(h) +
      Uniform(0.0001, 0.01).logPdf(sigma)
    case White(s) =>
      InverseGamma(3.0, 0.5).logPdf(s)
  }).sum

  val prior = for {
    v <- InverseGamma(3.0, 0.5)
    h <- InverseGamma(3.0, 4.0)
    sigma <- Uniform(0.0001, 0.05)
    beta <- Applicative[Rand].replicateA(3, Gaussian(0.0, 5.0))
  } yield
    GaussianProcess.Parameters(
      MeanParameters.plane(DenseVector(beta.toArray)),
      Vector(KernelParameters.se(h, sigma), KernelParameters.white(v))
    )

  val trainingData = data.
    filter(_.name != test_sensor).
    map(t => t.obs map (y => GaussianProcess.Data(Two(t.lon, t.lat), y))).
    flatten

  // get iterations from command line argument
  val nIters: Int = args.lift(0).map(_.toInt).getOrElse(100000)

  val iters = Mcmc.sample(trainingData, priorKernel,
    Gaussian(1.0, 1.0), proposal(0.05), dist, prior.draw).
    steps.
    take(nIters)

  val out = new java.io.File("data/temperature-mcmc.csv")
  val writer = out.asCsvWriter[List[Double]](rfc.withHeader("h", "sigma", "sigma_y", "beta_0", "beta_1", "beta_2"))

  def formatParams(p: GaussianProcess.Parameters): List[Double] = {
    p.kernelParameters.flatMap(_.toList).toList ::: p.meanParameters.toList
  }

  // write iters to file
  while (iters.hasNext) {
    writer.write(formatParams(iters.next))
  }

  writer.close()
}

object PredictTemperature extends App with TemperatureModel {
  val paramsFile = Paths.get("data/temperature-mcmc.csv")
  val paramsReader = paramsFile.asCsvReader[List[Double]](rfc.withHeader)
  val chain = paramsReader.
    collect {
      case Right(a) => GaussianProcess.Parameters(
        MeanParameters.plane(DenseVector(a(3), a(4), a(5))),
        Vector(KernelParameters.se(a(0), a(1)), KernelParameters.white(a(2)))
      )
    }.
    drop(1000).
    toStream.
    zipWithIndex.
    filter { case (_, i) => i % 20 == 0 }.
    map(_._1)

  // calculate mean of parameters
  val params: GaussianProcess.Parameters = chain.
    reduce { (a, b) =>
      GaussianProcess.Parameters(
        a.meanParameters add b.meanParameters,
        (a.kernelParameters zip b.kernelParameters).
          map { case (x, y) => x add y })
    }.
    map(p => p / chain.size)

  println(params)

  // read test locations
  val locationFile = Paths.get("data/spatial_data/test_location.csv")
  val locReader = locationFile.asCsvReader[(String, Double, Double)](rfc.withHeader)
  val locations = locReader.
    collect {
      case Right(a) => Two(a._2, a._3)
    }.
    toVector

  // perform prediction at test location for each day
  // TODO: Remove test locations
  val predictions = data.
    groupBy(_.date).
    map { case (d, ts) =>
      (d, ts.map( t => (GaussianProcess.Data(Two(t.lon, t.lat), t.
                                               obs.get))))
    }.
    toVector.
    flatMap { case (date, obs) =>
      val fitted = Predict.fit(locations, obs, dist, params)
      Predict.predict(fitted.map(_._2), 0.8).map(x => (date, x))
    }

  // write test prediction
  val out = new java.io.File("data/temperature_test_predictions_gp.csv")
  out.writeCsv(predictions, rfc.withHeader("day", "mean", "lower", "object"))
}
