import org.scalatest._
import prop._
// import org.scalactic.Equality
import org.scalacheck.Gen
import breeze.linalg.{DenseVector, DenseMatrix, diag}
import com.github.jonnylaw.gp._

// TODO: Write test for vector to params and param to vector

trait ParamGens {
  val denseVector: Int => Gen[DenseVector[Double]] = (n: Int) =>
    Gen
      .containerOfN[Array, Double](n, Gen.choose(-10.0, 10.0))
      .map(a => DenseVector(a))

  /**
    * Simulate a positive definite matrix with a given condition number
    */
  def symmetricPosDefMatrix(n: Int, c: Double) = {
    if (n > 2) {
      for {
        entries <- Gen.containerOfN[Array, Double](n - 2, Gen.choose(1.0, c))
        d = diag(DenseVector(Array(1.0, c) ++ entries))
        u <- denseVector(n)
        t = 2 / (u.t * u)
        i = DenseMatrix.eye[Double](n)
      } yield (i - t * u * u.t) * d * (i - t * u * u.t)
    } else {
      for {
        u <- denseVector(n)
        d = diag(DenseVector(1.0, c))
        t = 2 / (u.t * u)
        i = DenseMatrix.eye[Double](n)
      } yield (i - t * u * u.t) * d * (i - t * u * u.t)
    }
  }

  val se = for {
    sigma <- Gen.choose(2.0, 10.0)
    h <- Gen.choose(2.0, 10.0)
  } yield KernelParameters.se(h, sigma)

  val white = for {
    sigma <- Gen.choose(2.0, 10.0)
  } yield KernelParameters.white(sigma)

  val plane = for {
    beta <- denseVector(3)
  } yield MeanParameters.plane(beta)

  val params = for {
    kp <- Gen.containerOfN[Vector, KernelParameters](3, Gen.oneOf(se, white))
    mp <- plane
  } yield GaussianProcess.Parameters(mp, kp)
}

class Parameters
    extends PropSpec
    with GeneratorDrivenPropertyChecks
    with Matchers
    with ParamGens {

  property(
    "Parameters should convert to DenseVector and back given specification") {
    forAll(params) { (p: GaussianProcess.Parameters) =>
      assert(KernelParameters.arrayToParams(p, KernelParameters.paramsToArray(p)) === p)
    }
  }
}
