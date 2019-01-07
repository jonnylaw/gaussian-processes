package com.github.jonnylaw.gp

import com.cibo.evilplot.plot._
import com.cibo.evilplot.plot.aesthetics.DefaultTheme._
import com.cibo.evilplot.numeric.Point
import com.cibo.evilplot.plot.renderers._
import com.cibo.evilplot.geometry._
import breeze.stats.distributions.Gaussian
import cats.implicits._

object Plot {

  def scatterPlot(ys: Vector[GaussianProcess.Data]) = {
    val data = ys.map {
      case GaussianProcess.Data(One(x), y) =>
        Point(x, y)
      case _ => throw new Exception("Scatter plot can only display one dimensional data")
    }

    ScatterPlot(data)
  }

  /**
    * Plot the mean and intervals of a Gaussian Process
    */
  def gpPlot(gp: Vector[(Location[Double], Gaussian)]) = {
    val (mean, lower, upper) = gp.map {
      case (One(x), g) =>
        (Point(x, g.mean),
         Point(x, Summarise.getInterval(g.mean, g.variance, 0.05)),
         Point(x, Summarise.getInterval(g.mean, g.variance, 0.95)))
      case _ =>
        throw new Exception("GP plot can only display one dimensional data")
    }.unzip3

    Overlay(LinePlot(mean),
            LinePlot(lower,
                     pathRenderer = PathRenderer.default(
                       lineStyle = LineStyle.DashDot.some).some),
            LinePlot(upper,
                     pathRenderer = PathRenderer.default(
                       lineStyle = LineStyle.DashDot.some).some))
    .standard()
  }

  def ppPlot(gps: Vector[Vector[(Location[Double], Gaussian)]]) = {
    val toPlot = gps.map(gp => LinePlot(gp.map {
      case (One(x), g) => Point(x, g.mean)
      case _ => throw new Exception("PP plot can only display one dimensional data") }))

    Overlay(toPlot: _*).
      standard()
  }
}
