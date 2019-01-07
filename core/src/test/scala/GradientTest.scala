import org.scalatest._
import prop._
import breeze.linalg.{cholesky}
import org.scalactic.Equality
import com.github.jonnylaw.gp._

class Gradient
    extends PropSpec
    with GeneratorDrivenPropertyChecks
    with Matchers
    with ParamGens {

  val input = for {
    m <- symmetricPosDefMatrix(2, 100)
    d <- denseVector(2)
  } yield (m, d)

  ignore("forward solve should be equivalent to backslash in breeze") {
    forAll(input) {
      case (m, b) =>
        val l = cholesky(m)

        val x1 = Predict.forwardSolve(l, b)
        val x2 = m \ b

        assert(x1 == x2)
    }
  }
}
