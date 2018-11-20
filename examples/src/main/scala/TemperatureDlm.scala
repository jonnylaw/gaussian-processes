package examples

import gp.core._
import breeze.stats.distributions._
import breeze.stats.mean
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

// TODO: Read in parameter posterior distribution from chapter 4 temperature model
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

  val seasonalDlm = Dlm.polynomial(1) |+| Dlm.seasonal(24, 3)
  val mvDlm = List.fill(8)(seasonalDlm).reduce(_ |*| _)
  val model = DlmGp.Model(None, mvDlm, dist)

  // create a vector of time series grouped by location
  val trainingData: Vector[(Location[Double], Vector[Data])] = data.
    groupBy(a => (a.lon, a.lat)).
    map { case (l, ts) =>
      (Two(l._1, l._2), ts.groupBy(_.date).
         map { case (t, d) => Data(t.getMillis, DenseVector(d.map(_.obs).toArray)) }.toVector)
    }.
    toVector

  val priorV = InverseGamma(3.0, 3.0)
  val priorW = InverseGamma(3.0, 3.0)

  val prior = for {
    v <- Applicative[Rand].replicateA(8, priorV)
    w <- Applicative[Rand].replicateA(56, priorW)
    m0 <- Applicative[Rand].replicateA(56, Gaussian(0.0, 1.0))
    c0 <- Applicative[Rand].replicateA(56, InverseGamma(3.0, 3.0))
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
                              model, trainingData, prior.draw)

  def formatParameters(s: FitDlmGp.State): List[Double] = ???

  Streaming.writeParallelChain(
    iters, 2, 100000, "examples/data/temperature_dlmgp_share_state", formatParameters).
    runWith(Sink.onComplete { s =>
              println(s)
              system.terminate()
            })
}

// Fit a zero mean covariance function to the residuals of the Temperature thingy
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
    Vector(KernelParameters.se(0.5, 1.0), KernelParameters.white(1.0))
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
  val priorSigmaY = InverseGamma(3.0, 3.0)

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
