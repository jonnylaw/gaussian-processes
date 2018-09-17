package examples

import gp.core._
import breeze.stats.distributions._
import breeze.linalg.DenseVector
import breeze.numerics.lbeta
import core.dlm.model._
import akka.actor.ActorSystem
import akka.stream._
import scaladsl._
import math._
import com.cibo.evilplot.plot.aesthetics.DefaultTheme._
import Diagnostics._

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

  val params = DenseVector(0.8, 1.0, 0.2)
  val sims = arSims(params).steps.take(500).toVector

  /**
    * Transform the parameters to be constrained 
    */
  def constrain(p: DenseVector[Double]): Vector[Hmc.Parameter] = {
    val phi = Hmc.bounded(0.1, 0.9)(p(0))
    val mu = Hmc.unbounded(p(1))
    val sigma = Hmc.boundedBelow(0)(p(2))

    Vector(phi, mu, sigma)
  }

  // sample from the prior distributions
  val unconstrained = for {
    phi <- new Beta(5, 2)
    mu <- Gaussian(1.0, 1.0)
    sigma <- InverseGamma(10, 2)
  } yield DenseVector(Dglm.logit(phi), mu, log(sigma))


  /**
    * Calculate the log-prior in the transformed parameter space
    */
  def logPrior(p: DenseVector[Double]): Double = {
    val ps = constrain(p)
    val phi = ps(0).constrained
    val mu = ps(1).constrained
    val sigma = ps(2).constrained

    val priorPhi = new Beta(3, 4).logPdf(phi)
    val priorMu = Gaussian(0, 1).logPdf(mu)
    val priorSigma = new CauchyDistribution(0.0, 2.0).logPdf(sigma)

    priorPhi + priorMu + priorSigma + ps.map(x => abs(x.logJacobian)).sum
  }

  def ll(alphas: Vector[Double])(p: DenseVector[Double]) = {
    val ps = constrain(p)
    val phi = ps(0).constrained
    val mu = ps(1).constrained
    val sigma = ps(2).constrained

    val initsd = math.sqrt(math.pow(sigma, 2) / (1 - math.pow(phi, 2)))

    Gaussian(mu, initsd).logPdf(alphas.head) +
      (alphas.drop(2) zip alphas.tail.init).map {
        case (a1, a0) =>
          Gaussian(mu + phi * (a0 - mu), sigma).logPdf(a1)
      }.sum
  }

  /**
    * Gradient of the un-normalised log-posterior with
    * respect to the parameters, need to consider the derivative of the 
    * log-Jacobian 
    */
  def grad(alphas: Vector[Double])(
    p: DenseVector[Double]): DenseVector[Double] = {
    val ps = constrain(p)
    val phi = ps(0).constrained
    val mu = ps(1).constrained
    val sigma = ps(2).constrained

    val n = alphas.size

    val alphaPhi = 3.0
    val betaPhi = 4.0

    // cauchy
    val lSigma = 0.0
    val gamma = 2.0

    val muMu = 0.0
    val sigmaMu = 1.0

    val ssa = (alphas.init, alphas.tail).zipped.map {
      case (a0, a1) => (mu - a0) * (a1 - mu - phi * (a0 - mu))
    }.sum

    val dphi = (alphaPhi - 1.0)/ phi + (1.0 - betaPhi) / (1.0 - phi) + (1.0 / (sigma * sigma)) * ssa + ps(0).derivative

    val ssb = (alphas.init, alphas.tail).zipped.map {
      case (a0, a1) => a1 - (mu + phi * (a0 - mu))
    }.sum

    val dmu = -(1.0 / sigmaMu * sigmaMu) * (mu - muMu) + ((1.0 - phi) / pow(sigma, 2)) * ssb + ps(1).derivative

    val ssc = (alphas.init, alphas.tail).zipped.map {
      case (a0, a1) => math.pow(a1 - (mu + phi * (a0 - mu)), 2)
    }.sum

    val dsigma = -2.0/(Pi*math.pow(gamma + sigma - lSigma, 2)) - (n.toDouble / sigma) + (1.0 / pow(sigma, 3)) * ssc + ps(2).derivative

    DenseVector(dphi, dmu, dsigma)
  }
}

object Ar1Hmc extends App with Ar1Model {
  // implicit val system = ActorSystem("ar1_hmc")
  // implicit val materializer = ActorMaterializer()

  // implicit val basis = RandBasis.withSeed(103)

  val pos = (p: DenseVector[Double]) => ll(sims)(p) + logPrior(p)

  val iters = Hmc(3, 10, 0.0025, grad(sims), pos).
    sample(unconstrained.draw)

  def format(s: HmcState): Vector[Double] = {
    constrain(s.theta).map(_.constrained).toVector// ++ List(s.accepted.toDouble)
  }

  val chain = iters.
    steps.
    map(format).
    drop(10000).
    take(10000).
    toVector

  // (1 to 2).par.
  //   map { i =>
  //     val chain = iters.steps.map(format).take(10000).toVector

  diagnostics(chain).
    render().
    write(new java.io.File(s"figures/ar_hmc.png"))

  // Streaming
  //   .writeParallelChain(iters, 2, 1000, "examples/data/ar1_hmc", format)
  //   .runWith(Sink.onComplete(_ => system.terminate()))
}

object Ar1Mh extends App with Ar1Model {
  val pos = (p: DenseVector[Double]) => ll(sims)(p) + logPrior(p)
  val prop = (p: DenseVector[Double]) => {
    val z = DenseVector.rand(p.length, Gaussian(0, 0.1))
    Rand.always(p + z)
  }

  val chain = MarkovChain.metropolis(unconstrained.draw, prop)(ll(sims)).
    steps.
    take(10000).
    map(s => constrain(s).map(_.constrained)).
    toVector

  diagnostics(chain).
    render().
    write(new java.io.File(s"figures/ar_mh.png"))
}

object Ar1Da extends App with Ar1Model {
  implicit val system = ActorSystem("ar1_nuts")
  implicit val materializer = ActorMaterializer()

  implicit val basis = RandBasis.withSeed(4)

  val pos = (p: DenseVector[Double]) => logPrior(p) + ll(sims)(p)

  val m = DenseVector.ones[Double](3)
  val iters = DualAverage(10.0, 0.65, m, 1000, grad(sims), pos).
    sample(unconstrained.draw)

  def format(s: DualAverageState): List[Double] = {
    constrain(s.theta).map(_.constrained).toList ++
    List(s.accepted.toDouble, exp(s.logeps))
  }

  Streaming
    .writeParallelChain(iters, 2, 10000, "examples/data/ar1_da", format)
    .runWith(Sink.onComplete(_ => system.terminate()))
}

object Ar1Nuts extends App with Ar1Model {
  implicit val system = ActorSystem("ar1_nuts")
  implicit val materializer = ActorMaterializer()

  implicit val basis = RandBasis.withSeed(4321)

  val pos = (p: DenseVector[Double]) => logPrior(p) + ll(sims)(p)

  def format(s: DualAverageState): List[Double] = 
    constrain(s.theta).map(_.constrained).toList

  val m = DenseVector.ones[Double](3)
  val iters = Nuts(m, 0.65, 1000, grad(sims), pos).
    sample(unconstrained.draw)

  Streaming
    .writeParallelChain(iters, 2, 10000, "examples/data/ar1_nuts", format)
    .runWith(Sink.onComplete(_ => system.terminate()))
}
