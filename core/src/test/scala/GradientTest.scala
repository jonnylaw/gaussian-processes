import org.scalatest._
import prop._
import breeze.linalg.{cholesky}
import org.scalactic.Equality
import gp.core._

class Gradient
    extends PropSpec
    with GeneratorDrivenPropertyChecks
    with Matchers
    with ParamGens {

  val xs = GaussianProcess.samplePoints(-10.0, 10.0, 30).map(One.apply)
  val ys = (p: GaussianProcess.Parameters) =>
    GaussianProcess.draw(xs, dist, p).draw
  val data = (p: GaussianProcess.Parameters) =>
    (xs.map(_.x) zip ys(p).data.toVector).map {
      case (x, y) => GaussianProcess.Data(One(x), y)
  }
  val dist = Location.euclidean _

  // property(
  //   "simplification of gradient should be equivalent to naive implementation") {
  //   forAll(params) { (p: GaussianProcess.Parameters) =>
  //     GaussianProcess.mllGradient(data(p), dist)(p) == GaussianProcess
  //       .mllGradientNaive(data(p), dist)(p)
  //   }
  // }

  val input = for {
    m <- symmetricPosDefMatrix(2, 100)
    d <- denseVector(2)
  } yield (m, d)

  property("forward solve should be equivalent to backslash in breeze") {
    forAll(input) {
      case (m, b) =>
        val l = cholesky(m)

        val x1 = Predict.forwardSolve(l, b)
        val x2 = m \ b

        assert(x1 == x2)
    }
  }
}
