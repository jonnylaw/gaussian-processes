package examples

import gp.core._
import breeze.stats.distributions.{Gaussian, Uniform, Gamma, Rand}
import breeze.linalg.{diag, DenseVector}
import cats._
import cats.implicits._
import kantan.csv._
import kantan.csv.ops._
import java.nio.file.Paths
import core.dlm.model._
import akka.actor.ActorSystem
import akka.stream._
import scaladsl._

trait TestModel {
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

  val params = GaussianProcess.Parameters(
    MeanParameters.plane(DenseVector(-2.0, 0.05)),
    Vector(KernelParameters.se(1.0, 2.0), KernelParameters.white(0.09))
  )
  
  val covFn = KernelFunction.apply(params.kernelParameters)
  val dist = Location.euclidean _

  val prior = for {
    v <- InverseGamma(3.0, 0.5)
    h <- InverseGamma(3.0, 4.0)
    sigma <- Rand.always(2.0) // InverseGamma(0.001, 0.001)
    alpha <- Gaussian(0.0, 5.0)
    beta <- Gaussian(0.0, 5.0)
  } yield GaussianProcess.Parameters(
    MeanParameters.plane(DenseVector(alpha, beta)),
    Vector(KernelParameters.se(h, sigma), KernelParameters.white(v))
  )
}

object SimulateGp extends App with TestModel {
  val xs = GaussianProcess.samplePoints(-10.0, 10.0, 30).map(One.apply)
  val ys = GaussianProcess.draw(xs, dist, params).draw

  val out = new java.io.File("data/simulated_gp.csv")
  out.writeCsv(xs.map(_.x) zip ys.data.toVector, rfc.withHeader("x", "y"))
}

object FitGp extends App with TestModel {
  // read in data
  val rawData = Paths.get("data/simulated_gp.csv")
  val reader = rawData.asCsvReader[List[Double]](rfc.withHeader)
  val data = reader.
    collect { 
      case Right(a) => GaussianProcess.Data(One(a.head), a(1))
    }.
    toVector

  // create a regular grid of points to draw from the GP posterior
  implicit val integralD = scala.math.Numeric.DoubleAsIfIntegral
  val testPoints = Vector.range(-10.0, 10.0, 0.01).map(One(_))

  // draw from the GP posterior at values of testPoints
  val fitted = Predict.fit(testPoints ++ data.map(_.x), data, dist, params)

  // combine 
  val out = (testPoints.map(_.x) zip fitted.map(_.mu) zip
    Summarise.getIntervals(DenseVector(fitted.map(_.mu).toArray),
      diag(DenseVector(fitted.map(_.sigma).toArray)), 0.975)).
    map { case ((a, b), (c, d)) => (a, b, c, d) }

  val outFile = new java.io.File("data/fitted_gp.csv")
  outFile.writeCsv(out, rfc.withHeader("x", "mean", "upper", "lower"))
}

object ParametersSimulatedGp extends App with TestModel {
  implicit val system = ActorSystem("fit_simulated_gp")
  implicit val materializer = ActorMaterializer()

  val rawData = Paths.get("data/simulated_gp.csv")
  val reader = rawData.asCsvReader[List[Double]](rfc.withHeader)
  val data = reader.
    collect { 
      case Right(a) => GaussianProcess.Data(One(a.head), a(1))
    }.
    toVector

  def proposal(delta: Double)(ps: Vector[KernelParameters]) = ps traverse {
    p => p match {
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
    }}

  def priorKernel(ps: Vector[KernelParameters]) = ps.map(p => p match {
    case SquaredExp(h, sigma) =>
      InverseGamma(0.001, 0.001).logPdf(h) +
      InverseGamma(0.001, 0.001).logPdf(sigma)
    case White(s) =>
      InverseGamma(0.001, 0.001).logPdf(s)
  }).sum

  // get iterations from command line argument
  val nIters: Int = args.lift(0).map(_.toInt).getOrElse(100000)

  val iters = Mcmc.sampleWithState(data, Gamma(3, 0.5),
    priorKernel, Gaussian(1.0, 1.0), proposal(0.05), dist, prior.draw)

  def format(s: Mcmc.State): List[Double] = 
    s.p.kernelParameters.flatMap(_.toList).toList ::: s.p.meanParameters.toList
  
  // write iters to file
  Streaming.
    writeParallelChain(iters, 2, 100000, "data/gpmcmc", format).
    runWith(Sink.onComplete(_ => system.terminate()))
}
