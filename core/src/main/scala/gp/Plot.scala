package com.github.jonnylaw.gp

import com.cibo.evilplot.plot._
import com.cibo.evilplot.plot.aesthetics.DefaultTheme._
import com.cibo.evilplot.numeric._
import com.cibo.evilplot.plot.renderers._
import com.cibo.evilplot.colors._
import com.cibo.evilplot.geometry._
import breeze.stats.distributions.Gaussian
import cats.implicits._

object Plot {

  /**
    * Plot a Gaussian process with a one-dimensional location, x
    * @param ys a vector containing GP Data, a product type with locations
    * and data
    * @return a scatter plot
    */
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
    * @param gp 
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
    val toPlot = gps.map { gp =>
      LinePlot(gp.map {
                 case (One(x), g) => Point(x, g.mean)
                 case _ =>
                   throw new Exception("PP plot can only display one dimensional location data") },
               pathRenderer = PathRenderer.default(color = HTMLNamedColors.red.copy(opacity = 0.5).some).some) }

    Overlay(toPlot: _*).
      standard()
  }

  /**
    * Plot a contour plot of a GP with two-dimensional location data
    */
  // def contourPlot(ys: Vector[GaussianProcess.Data]) = {
  //   val locationData = ys.map {
  //     case GaussianProcess.Data(Two(x1, x2), _) => Point(x1, x2)
  //     case _ =>    throw new Exception("surface plot can only display two dimensional location data")
  //   }
  //   val data = ys.map {
  //     case GaussianProcess.Data(Two(x1, x2), y) => Point3(x1, x2, y)
  //     case _ =>    throw new Exception("surface plot can only display two dimensional location data")
  //   }

  //   ContourPlot(locationData,
  //               surfaceRenderer = SurfaceRenderer.densityColorContours(data))
  // }
}
