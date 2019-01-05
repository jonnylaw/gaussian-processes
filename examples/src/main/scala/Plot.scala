package com.github.jonnylaw.gp

import com.cibo.evilplot.plot._
import com.cibo.evilplot.plot.aesthetics.DefaultTheme._
import com.cibo.evilplot.numeric.Point

object Plot {

  def scatterPlot(ys: Vector[Data]): Unit = {
    val data = ys.map {
      case (One(x), y) =>
        Point(x, y)
      case _ => throw new Exception("Scatter plot can only display one dimensional data")
    }

    ScatterPlot(data)
      .xAxis()
      .yAxis()
      .frame()
      .xLabel("x")
      .yLabel("y")
      .render()
  }
}
