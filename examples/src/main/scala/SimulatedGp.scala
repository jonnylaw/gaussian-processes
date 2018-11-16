package examples

import gp.core._
import breeze.stats.distributions.{Gaussian, Uniform, Gamma, Rand}
import breeze.linalg.{diag, DenseVector}
import cats._
import cats.implicits._
import kantan.csv._
import kantan.csv.ops._
import java.nio.file.Paths
import dlm.core.model._
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
    // MeanParameters.plane(DenseVector(-2.0, 0.5)),
    MeanParameters.zero,
    Vector(KernelParameters.se(3.0, 5.5), KernelParameters.white(0.5))
  )

  val covFn = KernelFunction.apply(params.kernelParameters)
  val dist = Location.euclidean _

  val prior = for {
    v <- InverseGamma(3.0, 0.5)
    h <- InverseGamma(3.0, 4.0)
    sigma <- Rand.always(2.0) // InverseGamma(0.001, 0.001)
    // alpha <- Gaussian(0.0, 5.0)
    // beta <- Gaussian(0.0, 5.0)
  } yield GaussianProcess.Parameters(
    // MeanParameters.plane(DenseVector(alpha, beta)),
    MeanParameters.zero,
    Vector(KernelParameters.se(h, sigma), KernelParameters.white(v))
  )
}

object SimulateGp extends App with TestModel {
  val xs = GaussianProcess.samplePoints(-10.0, 10.0, 300).map(One.apply)
  val ys = GaussianProcess.draw(xs, dist, params).draw.data.toVector

  val out = new java.io.File("data/simulated_gp.csv")
  out.writeCsv(xs.map(_.x) zip ys, rfc.withHeader("x", "y"))
}

// simulate p independent realisations of a Gaussian process from the test model
object SimulateGpReplicate extends App with TestModel {
  val xs = GaussianProcess.samplePoints(-10.0, 10.0, 30).map(One.apply)
  val p = 50

  val res: Vector[(Int, Double, Double)] = Vector.range(1, p).
    flatMap(_ => GaussianProcess.draw(xs, dist, params).draw.data.toVector).
    zip((1 to p).flatMap(i => Vector.fill(xs.size)(i))).
    zip(Vector.fill(p)(xs.map(_.x)).flatten).
    map { case ((y, i), x) => (i, x, y) }

  val out = new java.io.File("data/simulated_gp_replicate.csv")
  out.writeCsv(res, rfc.withHeader("replicate", "x", "y"))
}

object FitGp extends App with TestModel {
  // read in data
  val rawData = Paths.get("data/simulated_gp.csv")
  val reader = rawData.asCsvReader[List[Double]](rfc.withHeader)
  val data = reader.
    collect {
      case Right(a) => GaussianProcess.Data(One(a.head), a(1))
    }.
    toVector.
    zipWithIndex.
    filter { case (_, i) => i % 15 == 0 }.
    map(_._1)

  // create a regular grid of points to draw from the GP posterior
  implicit val integralD = scala.math.Numeric.DoubleAsIfIntegral
  val testPoints = Vector.range(-10.0, 10.0, 0.01).map(One(_))

  // draw from the GP posterior at values of testPoints
  val fitted = Predict.fit(testPoints ++ data.map(_.x), data, dist, params)

  // combine
  val out = (fitted.map { case (x, g) => (x.toVector.head, g.mean) } zip
    Summarise.getIntervals(DenseVector(fitted.map(_._2.mean).toArray),
      diag(DenseVector(fitted.map(_._2.variance).toArray)), 0.975)).
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
    toVector.
    zipWithIndex.
    filter { case (_, i) => i % 15 == 0 }.
    map(_._1)

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
      InverseGamma(3, 6).logPdf(h) +
      InverseGamma(3, 6).logPdf(sigma)
    case White(s) =>
      InverseGamma(10, 6).logPdf(s)
  }).sum

  // get iterations from command line argument
  val nIters: Int = args.lift(0).map(_.toInt).getOrElse(100000)

  val iters = Mcmc.sample(data, priorKernel,
                          Gaussian(1.0, 1.0), proposal(0.05), dist, prior.draw)

  def format(p: GaussianProcess.Parameters): List[Double] =
    p.kernelParameters.flatMap(_.toList).toList ::: p.meanParameters.toList

  // write iters to file
  Streaming.
    writeParallelChain(iters, 2, nIters, "data/gpmcmc", format).
    runWith(Sink.onComplete(_ => system.terminate()))
}

object PosteriorPredictive extends App with TestModel {
  val rawData = Paths.get("data/simulated_gp.csv")
  val reader = rawData.asCsvReader[List[Double]](rfc.withHeader)
  val data = reader.
    collect {
      case Right(a) => GaussianProcess.Data(One(a.head), a(1))
    }.
    toVector.
    zipWithIndex.
    filter { case (_, i) => i % 15 == 0 }.
    map(_._1)

  val rawChain = Paths.get("data/gpmcmc_0.csv")
  val paramReader = rawChain.asCsvReader[List[Double]](rfc.withHeader)
  val listParams = paramReader.
    collect {
      case Right(a) => GaussianProcess.Parameters(
        MeanParameters.zero,
        Vector(KernelParameters.se(a.head, a(1)), KernelParameters.white(a(2)))
      )
    }.
    toVector.
    drop(10000) // remove burn-in

  def discreteUniform(min: Int, max: Int) = new Rand[Int] {
    def draw = scala.util.Random.nextInt(max - min) + min
  }

  val p = 100

  // sample from the parameters and perform a fit
  val indices = discreteUniform(0, listParams.size).sample(p).toVector

  // create a regular grid of points to draw from the GP posterior
  implicit val integralD = scala.math.Numeric.DoubleAsIfIntegral
  val testPoints = Vector.range(-10.0, 10.0, 0.01).map(One(_))

  def predict(p: GaussianProcess.Parameters): Vector[List[Double]]  = {
    Predict.fit(testPoints, data, dist, p).
      map { case (x, g) => x.toVector.head :: g.mean :: p.kernelParameters.flatMap(_.toList).toList }
  }

  val parameters = indices.toVector.map(i => listParams(i))
  val res: Vector[List[Double]] = parameters.
    map(predict).
    zipWithIndex.
    flatMap { case (d, i) => d.map(l => i.toDouble :: l) }

  val out = new java.io.File("data/posterior_predictive.csv")
  out.writeCsv(res, rfc.withHeader("replicate", "x", "mean", "h", "sigma", "sigma_y"))
}
