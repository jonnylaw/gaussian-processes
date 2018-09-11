package examples

import gp.core._
import breeze.stats.distributions._
import breeze.linalg.DenseVector
import breeze.numerics.lbeta
import core.dlm.model._
import akka.actor.ActorSystem
import akka.stream._
import scaladsl._

trait Ar1Model {
  def arStep(p: DenseVector[Double])(a0: Double) = {
    val phi = p(0)
    val mu = p(1)
    val sigma = p(2)

    Gaussian(mu + phi * (a0 - mu), sigma)
  }

  def arSims(p: DenseVector[Double]): Process[Double] = {
    val phi = p(0)
    val mu = p(1)
    val sigma = p(2)

    val init =
      Gaussian(mu, math.sqrt(sigma * sigma / (1 - phi * phi))).draw
    MarkovChain(init)(arStep(p))
  }

  val params = DenseVector(0.8, 1.0, 0.3)
  val sims = arSims(params).steps.take(500).toVector
}

object ArHmc extends App with Ar1Model {
  implicit val system = ActorSystem("fit_simulated_gp")
  implicit val materializer = ActorMaterializer()

  def logPrior(p: DenseVector[Double]): Double = {
    val phi = p(0)
    val mu = p(1)
    val sigma = p(2)

    val priorSigma =
      if (sigma < 0) -1e99
      else new CauchyDistribution(1.0, 1.0).logPdf(sigma)

    val priorPhi =
      if (phi > 1 || phi < 0) -1e99
      else new Beta(3, 4).logPdf(phi)

    priorPhi +
      Gaussian(0, 1).logPdf(mu) +
      priorSigma
  }

  def ll(alphas: Vector[Double])(p: DenseVector[Double]) = {
    val phi = p(0)
    val mu = p(1)
    val sigma = p(2)

    val initsd = math.sqrt(math.pow(sigma, 2) / (1 - math.pow(phi, 2)))
    Gaussian(mu, initsd).logPdf(alphas.head) +
      (alphas.tail zip alphas.init).map {
        case (a1, a0) =>
          Gaussian(mu + phi * (a0 - mu), sigma).logPdf(a0)
      }.sum
  }

  // gradient of the log-posterior with respect to the parameters
  def grad(alphas: Vector[Double])(
      p: DenseVector[Double]): DenseVector[Double] = {
    val phi = p(0)
    val mu = p(1)
    val sigma = p(2)

    val n = alphas.size

    val alphaPhi = 3.0
    val betaPhi = 4.0

    // val l_sigma = 1.0

    val ssa = (alphas.init, alphas.tail).zipped.map {
      case (a0, a1) => (a0 - mu) * (a1 - (mu + phi * (a0 - mu)))
    }.sum

    val dphi = (alphaPhi - 1) * (betaPhi - 1) * (1 - 2 * phi) / (phi * (1 - phi) * lbeta(
      alphaPhi,
      betaPhi)) + 1 / (sigma * sigma) * ssa

    val ssb = (alphas.init, alphas.tail).zipped.map {
      case (a0, a1) => a1 - (mu + phi * (a0 - mu))
    }.sum

    val dmu = -mu - (phi - 1) / math.pow(sigma, 2) * ssb

    val ssc = (alphas.init, alphas.tail).zipped.map {
      case (a0, a1) => math.pow(a1 - (mu + phi * (a0 - mu)), 2)
    }.sum

    val dsigma = -n / sigma + (1 / math.pow(sigma, 3)) * ssc

    // Cauchy Prior derivative
    // 2*sigma/(sigma-l_sigma)

    DenseVector(dphi, dmu, dsigma)
  }

  implicit val basis = RandBasis.withSeed(2)

  val m = DenseVector.ones[Double](3)
  val pos = (p: DenseVector[Double]) => logPrior(p) + ll(sims)(p)
  val iters = Hmc(0.01, 0.65, m, 1000, grad(sims), pos).sample(params)

  def format(s: HmcState): List[Double] = {
    s.theta.data.toList ++ List(s.accepted.toDouble)
  }

  Streaming
    .writeParallelChain(iters, 2, 10000, "examples/data/ar1_hmc", format)
    .runWith(Sink.onComplete(_ => system.terminate()))
}

object Ar1Nuts extends App {}
