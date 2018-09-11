package examples

import gp.core._
import com.stripe.rainier.compute._
import com.stripe.rainier.core._
import com.stripe.rainier.sampler._
import kantan.csv._
import kantan.csv.ops._

/**
  * Use automatic differentiation and HMC for inference
  * for the AR(1) parameter posterior distribution
  */
object RanierAr extends App with Ar1Model {
  implicit val rng = ScalaRNG(4) // set a seed

  case class Parameters(phi: Real, mu: Real, sigma: Real)

  val prior = for {
    phi <- Beta(2, 5).param
    mu <- Normal(0, 1).param
    sigma <- Cauchy(0, 2).param
  } yield Parameters(phi, mu, sigma)

  def addTimePoint(params: RandomVariable[Parameters],
                   alphas: (Double, Double)): RandomVariable[Parameters] =
    for {
      p <- params
      _ <- Normal(p.mu + p.phi * (alphas._2 - p.mu), p.sigma).fit(alphas._1)
    } yield p

  val fullModel = sims.tail.zip(sims.init).foldLeft(prior)(addTimePoint)

  val model = for {
    p <- fullModel
  } yield Map("mu" -> p.mu, "phi" -> p.phi, "sigma" -> p.sigma)

  val iters = model.sample(HMC(10), 1000, 10000, 1)

  val out = new java.io.File("examples/data/ar1_rainier.csv")
  val headers = rfc.withHeader("mu", "phi", "sigma")

  out.writeCsv(iters.map(_.values), headers)
}
