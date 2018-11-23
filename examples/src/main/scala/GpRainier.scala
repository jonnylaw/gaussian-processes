package examples

import gp.core._
import com.stripe.rainier.compute._
import com.stripe.rainier.core._
import com.stripe.rainier.sampler._
import java.nio.file.Paths
import kantan.csv._
import kantan.csv.ops._

object FitSimulatedGpHmc extends App {
  val rawData = Paths.get("data/simulated_gp.csv")
  val reader = rawData.asCsvReader[List[Double]](rfc.withHeader)
  val data = reader.collect {
    case Right(a) => GaussianProcess.Data(One(a.head), a(1))
  }.toVector

  case class Parameters(h: Real,
                        sigma: Real,
                        sigmaY: Real,
                        alpha: Real,
                        beta: Real)

  val prior = for {
    v <- Gamma(0.001, 0.001).param
    h <- Gamma(0.001, 0.001).param
    sigma <- Gamma(0.001, 0.001).param
    alpha <- Normal(0.0, 5.0).param
    beta <- Normal(0.0, 5.0).param
  } yield Parameters(h, sigma, v, alpha, beta)

  // TODO: Write Rainier model for GP
  // def dist =

  // def model(ys: Vector[Data]) = for {
  //   _ = Normal(    h * exp(- (dist * dist) / (sigma * sigma) )

  // } yield p
}
