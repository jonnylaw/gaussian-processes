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

// TODO: Consider sampling on the unbounded space more carefully
// transformations
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

  /**
    * Transform the parameters to be constrained 
    */
  def constrain(p: DenseVector[Double]): Vector[Hmc.Parameter] = {
    val phi = Hmc.bounded(0, 1)(p(0))
    val mu = Hmc.unbounded(p(1))
    val sigma = Hmc.boundedBelow(0)(p(2))

    Vector(phi, mu, sigma)
  }

  /**
    * Calculate the log-prior in the transformed parameter space
    */
  def logPrior(p: DenseVector[Double]): Double = {
    val ps = constrain(p)
    println(s"Evaluting the prior with these parameters ${ps.map(_.constrained)}")

    val priorPhi = new Beta(3, 4).logPdf(ps(0).constrained) + ps(0).logJacobian
    val priorMu = Gaussian(0, 1).logPdf(ps(1).constrained) + ps(1).logJacobian
    val priorSigma = new CauchyDistribution(1.0, 1.0).logPdf(ps(2).constrained) + ps(2).logJacobian 

    priorPhi + priorMu + priorSigma
  }

  def ll(alphas: Vector[Double])(p: DenseVector[Double]) = {
    val ps = constrain(p)
    val phi = ps(0).constrained
    val mu = ps(1).constrained
    val sigma = ps(2).constrained

    val initsd = math.sqrt(math.pow(sigma, 2) / (1 - math.pow(phi, 2)))
    Gaussian(mu, initsd).logPdf(alphas.head) +
      (alphas.tail zip alphas.init).map {
        case (a1, a0) =>
          Gaussian(mu + phi * (a0 - mu), sigma).logPdf(a0)
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
    val lSigma = 1.0

    val ssa = (alphas.init, alphas.tail).zipped.map {
      case (a0, a1) => (a0 - mu) * (a1 - (mu + phi * (a0 - mu)))
    }.sum

    val dphi = (alphaPhi - 1) * (betaPhi - 1) * (1 - 2 * phi) / (phi * (1 - phi) * lbeta(alphaPhi, betaPhi)) + (1 / (sigma * sigma) * ssa) + ps(0).derivative

    val ssb = (alphas.init, alphas.tail).zipped.map {
      case (a0, a1) => a1 - (mu + phi * (a0 - mu))
    }.sum

    val dmu = -mu + ((1 - phi) / math.pow(sigma, 2)) * ssb + ps(1).derivative

    val ssc = (alphas.init, alphas.tail).zipped.map {
      case (a0, a1) => math.pow(a1 - (mu + phi * (a0 - mu)), 2)
    }.sum

    // Cauchy Prior derivative
    val dsigma = 2*sigma/(sigma-lSigma) -n / sigma + (1 / pow(sigma, 3)) * ssc + ps(2).derivative

    DenseVector(dphi, dmu, dsigma)
  }
}

object ArHmc extends App with Ar1Model {
  // implicit val system = ActorSystem("fit_simulated_gp")
  // implicit val materializer = ActorMaterializer()

  implicit val basis = RandBasis.withSeed(2)

  val pos = (p: DenseVector[Double]) => logPrior(p) + ll(sims)(p)

  val unconstrained =
    DenseVector(Dglm.logit(params(0)), params(1), log(params(2)))

  println(s"initial parameters $unconstrained")
  println(s"constrained parameters: ${constrain(unconstrained)}")

  val iters = Hmc(3, 5, 0.05, grad(sims), pos).
    sample(unconstrained)

  def format(s: HmcState): List[Double] = {
    constrain(s.theta).map(_.constrained).toList ++ List(s.accepted.toDouble)
  }

  iters.
    steps.
    take(100).
    map(s => constrain(s.theta).map(_.constrained)).
    foreach(println)

  // Streaming
  //   .writeParallelChain(iters, 2, 1000, "examples/data/ar1_hmc", format)
  //   .runWith(Sink.onComplete(_ => system.terminate()))
}

object Ar1Nuts extends App with Ar1Model {
  implicit val basis = RandBasis.withSeed(4)

  val pos = (p: DenseVector[Double]) => logPrior(p) + ll(sims)(p)

  val unconstrained =
    DenseVector(Dglm.logit(params(0)), params(1), log(params(2)))

  val iters = Nuts(3, 0.05, 1000, grad(sims), pos, 0.5).
    sample(unconstrained)

  def format(s: HmcState): List[Double] = {
    s.theta.data.toList ++ List(s.accepted.toDouble)
  }

  iters.
    steps.
    take(100).
    map(s => constrain(s.theta).map(_.constrained)).
    foreach(println)
}

