package examples

import gp.core._
import com.stripe.rainier.compute._
import com.stripe.rainier.core._
import com.stripe.rainier.sampler._
import java.nio.file.Paths
import kantan.csv._
import kantan.csv.ops._

object HmcParameters extends App with TestModel {
  val rawData = Paths.get("examples/data/simulated_gp.csv")
  val reader = rawData.asCsvReader[List[Double]](rfc.withHeader)
  val data = reader.collect {
    case Right(a) => (a.head, a(1))
  }.toVector

  def format(s: HmcState): List[Double] = 
    s.theta.data.toList ++ List(s.accepted.toDouble)

  case class Parameters(h: Double, sigma: Double, sigmaY: Double)

  // def priorKernel(ps: Vector[KernelParameters]) =
  //   for {
  //     h <- Gamma(0.001, 0.001).param
  //     sigma <- Gamma(0.001, 0.001).param
  //     sigmaY <- Gamma(0.001, 0.001).param
  //   } yield Parameters(h, sigma, sigmaY)

  // def model(ys: Vector[(Double, Double)], p: GaussianProcess.Parameters) = {
  //   val kxx = KernelFunction.buildCov(
  //     xs,
  //     KernelFunction.apply(p.kernelParameters),
  //     dist)

  //   Normal(
  // }
}
