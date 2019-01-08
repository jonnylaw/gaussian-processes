---
layout: docs
title: "Parameter Inference"
---

# Hamiltonian Monte Carlo

Hamiltonian Monte Carlo (HMC) methods can be used to determine the posterior
distribution of the hyper-parameters. HMC uses the gradient of the log-posterior
in order to make more efficient proposals. 

<!-- ```tut:silent -->
<!-- import breeze.linalg.DenseMatrix -->

<!-- val alpha = 2.0 -->
<!-- val beta = 2.0 -->
<!-- val priorSigma = GradDist.gamma(Gamma(alpha, 1.0 / beta)) -->
<!-- val priorH = GradDist.gamma(Gamma(alpha, 1.0 / beta)) -->
<!-- val priorSigmaY = GradDist.gamma(Gamma(alpha, 1.0 / beta)) -->
<!-- val prior = Vector(priorSigma, priorH, priorSigmaY) -->

<!-- KernelParameters.sampleHmc(observed, dist, init, -->
<!--   prior, DenseMatrix.eye[Double](3), 5, 0.2). -->
<!--   steps. -->
<!--   take(10). -->
<!--   toVector. -->
<!--   foreach(println) -->
<!-- ``` -->

The parameter diagnostics for the HMC method with 5 leapfrog steps and step-size 0.2

<!-- ```tut:silent -->
<!-- Diagnostics.diagnostics(hmcIters.map { ps => ps.toMap }). -->
<!--   render(). -->
<!--   write(new java.io.File("docs/src/main/resources/figures/parameters_weakly_informative_gp.png")) -->
<!-- ``` -->
