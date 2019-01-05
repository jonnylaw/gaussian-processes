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
      .xAxis()
      .yAxis()
      .frame()
      .xLabel("x")
      .yLabel("y")
  }

  /**
    * Plot the mean and intervals of a Gaussian Process
    */
  def gpPlot(gp: Vector[(Location[Double], Gaussian)]) = {
    val (mean, lower, upper) = gp.map {
      case (One(x), g) =>
        (Point(x, g.mean), Point(x, Summarise.getInterval(g.mean, g.variance, 0.05)),
         Point(x, Summarise.getInterval(g.mean, g.variance, 0.95)))
      case _ =>
        throw new Exception("GP plot can only display one dimensional data")
    }.unzip3

    LinePlot(mean)
      .xAxis()
      .yAxis()
      .frame()
      .xLabel("x")
      .yLabel("y")

    Overlay(LinePlot(lower,
                     pathRenderer = PathRenderer.default(
                       lineStyle = LineStyle.DashDot.some).some)
            .xAxis()
            .yAxis()
            .frame()
            .xLabel("x")
            .yLabel("y"))

    Overlay(LinePlot(upper,
                     pathRenderer = PathRenderer.default(
                       lineStyle = LineStyle.DashDot.some).some)
            .xAxis()
            .yAxis()
            .frame()
            .xLabel("x")
            .yLabel("y"))
  }
}
