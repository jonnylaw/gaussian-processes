package gp.core

import GaussianProcess._
import breeze.stats.distributions._
import breeze.linalg.{DenseVector, DenseMatrix, cholesky, diag, inv}
import cats.implicits._
import com.github.fommil.netlib.BLAS.{getInstance => blas}
import cats.data.Kleisli

/**
  * Use Gibbs Sampling and Metropolis-Hastings
  * to determine the state and parameters of a model
  * group by locations in data
  */
object Mcmc {

  /**
    * The state of the MCMC sampler, containing the currently sampled
    * values of the parameters and values of the function
    * @param sigma the measurement noise
    * @param p the current parameters
    * @param fx realisations of the function we are learning
    */
  case class State(
    sigma: Double,
    p:     Parameters,
    fx:    Vector[Data])

  /**
    * Join two vectors with a condition
    * @param xs a vector containing a type A
    * @param ys a vector containing another type b
    * @param cond the condition 
    */
  def innerJoin[A, B](
    xs:   Vector[A],
    ys:   Vector[B],
    cond: (A, B) => Boolean): Vector[(A, B)] = {

    for {
      x <- xs
      y <- ys
      if cond(x, y)
    } yield (x, y)
  }

  /**
    * Sample the observation variance from a Gamma distributoin
    * Gibbs step using a conditionally-conjugate distribution
    * @param prior the conditional conjugate prior distribtion for the 
    * measurement noise variance
    * @param ys a vector of data points
    * @param fx the currently sampled value of the function state
    * @return a distribution over the measurement noise
    */
  def samplePrecY(
    prior: Gamma,
    ys:    Vector[Data],
    fx:    Vector[Data]) = {

    val ssy = innerJoin(ys, fx,
      (a: Data, b: Data) => a.x === b.x).map {
      case (y, f) => (y.y - f.y) * (y.y - f.y)
    }.sum

    val shape = prior.shape + ys.size * 0.5
    val scale = prior.scale + 0.5 * ssy

    Gamma(shape, scale)
  }

 /**
   * Backsolve an upper-triangular linear system
   * with a single RHS
   *
   * @param A An upper-triangular matrix
   * @param y A single vector RHS
   *
   * @return The solution, x, of the linear system A x = y
   */
  def backSolve(
    A: DenseMatrix[Double],
    y: DenseVector[Double]): DenseVector[Double] = {
    val yc = y.copy
    blas.dtrsv("U", "N", "N", A.cols, A.toArray,
      A.rows, yc.data, 1)
    yc
  }

  /** 
    * Backsolve an upper-triangular linear system
    * with multiple RHSs
    * 
    * @param A An upper-triangular matrix
    * @param Y A matrix with columns corresponding to RHSs
    * 
    * @return Matrix of solutions, X, to the linear system A X = Y
    */
  def backSolve(
    A: DenseMatrix[Double],
    Y: DenseMatrix[Double]): DenseMatrix[Double] = {
    val yc = Y.copy
    blas.dtrsm("L", "U", "N", "N", yc.rows, yc.cols, 1.0, A.toArray, A.rows, yc.data, yc.rows)
    yc
  }

  /**
    * Sample the Gaussian Process at the existing locations
    * using a Gibbs step (is this really necessary?)
    * @param ys a vector of data
    * @param xs a vector of distinct locations
    * @param dist a distance function
    * @param sigma observation variance
    * @return a function from State => Rand[State]
    */
  def sampleState(
    ys:    Vector[Data],
    xs:    Vector[Location[Double]],
    dist:  (Location[Double], Location[Double]) => Double,
    sigma: Double,
    p:     Parameters
  ) = {
    
    // val sample: Vector[Rand[Double]] = Predict.fit(xs, ys, dist, s.p)
    // val newfx = (sample.sequence.draw zip xs) map { case (f, x) => Data(x, f) }

    // Rand.always(s.copy(fx = newfx))

    val covFn = KernelFunction.apply(p.kernelParameters)
    val meanFn = MeanFunction.apply(p.meanParameters)

    val nugget = diag(DenseVector.fill(xs.size)(1e-6))
    val kxx = KernelFunction.buildCov(xs, covFn, dist) + nugget

    val invKx = inv(kxx)
    val obsPrec = diag(DenseVector.fill(xs.size)(1.0 / sigma))
    val prec = invKx + obsPrec

    val y = DenseVector(ys.map(_.y).toArray)
    
    val mean = prec \ (invKx * DenseVector(xs.map(meanFn).toArray) + obsPrec * y)

    val z = DenseVector.fill(mean.size)(Gaussian(0, 1).draw)
    val fxSample = mean + cholesky(prec) \ z
    val newfx = (fxSample.data.toVector zip xs) map { case (f, x) => Data(x, f) }

    Rand.always(newfx)
  }

  /**
    * Sample the Gaussian Process from the prior, this will result in slow mixing
    */
  def sampleStatePrior(
    xs:   Vector[Location[Double]],
    dist: (Location[Double], Location[Double]) => Double) = { s: State =>

    GaussianProcess.drawData(xs, dist, s.p) map (newfx => s.copy(fx = newfx))
  }

  /**
    * 
    */
  def step(
    ys:         Vector[Data],
    xs:         Vector[Location[Double]],
    priorSigma: Gamma,
    priorKern:  Vector[KernelParameters] => Double,
    priorMean:  Gaussian,
    propKern:   Vector[KernelParameters] => Rand[Vector[KernelParameters]],
    dist:       (Location[Double], Location[Double]) => Double) = { st: State =>

    for {
      fx <- sampleState(ys, xs, dist, st.sigma, st.p)
      precY <- samplePrecY(priorSigma, ys, fx)
      kp <- KernelParameters.sample(ys, dist, priorKern, propKern)(st.p)
      mp <- MeanParameters.sample(ys, xs, dist, priorMean)(st.p)
    } yield State(math.sqrt(1.0 / precY), Parameters(mp.meanParameters, kp.kernelParameters), fx)
  }

  /**
    * Perform Gibbs sampling for the Gaussian Process
    * @param ys a collection of multivariate observations of a spatial process
    * @param priorSigma the conditional conjugate prior distribution of 
    * the measurement noise
    * @param prior the prior distribution of the covariance function parameters
    * @param priorMeanI
    * @param proposal the (symmetric) proposal distribution for the covariance function
    * parameters
    * @param initP the initial value of the parameters of the Gaussian Process
    */
  def sampleWithState(
    ys:         Vector[Data],
    priorSigma: Gamma,
    priorKern:  Vector[KernelParameters] => Double,
    priorMean:  Gaussian,
    propKern:   Vector[KernelParameters] => Rand[Vector[KernelParameters]],
    dist:       (Location[Double], Location[Double]) => Double,
    initP:      Parameters
  ): Process[State] = {

    val xs = ys.map(_.x).distinct
    val init = for {
      initState <- GaussianProcess.draw(xs, dist, initP)
    } yield State(priorSigma.draw, initP, Predict.buildData(xs, initState.data.toVector))

    MarkovChain(init.draw)(step(ys, xs, priorSigma, priorKern,
      priorMean, propKern, dist))
  }

  /**
    * Sample the Parameters 
    */
  def sample(
    ys:         Vector[Data],
    priorKern:  Vector[KernelParameters] => Double,
    priorMean:  Gaussian,
    propKern:   Vector[KernelParameters] => Rand[Vector[KernelParameters]],
    dist:       (Location[Double], Location[Double]) => Double,
    initP:      Parameters
  ) = {

    val xs = ys.map(_.x).distinct

    val step = Kleisli(KernelParameters.sample(ys, dist, priorKern, propKern)) compose
      Kleisli(MeanParameters.sample(ys, xs, dist, priorMean))

    MarkovChain(initP)(step.run)
  }
}
