package examples

import java.io.{FileInputStream, DataInputStream}
import java.nio.file.Paths
import java.util.zip.GZIPInputStream
import breeze.linalg._

/**
  * Use Gaussian Process Classification to determine the labels of handwritten 
  * digits from the MNIST dataset
  */
trait ReadMnist {
  val trainingImagesFile = Paths.get("data/mnist/", "train-images-idx3-ubyte.gz")
  val trainingImages = new DataInputStream(
    new GZIPInputStream(
      new FileInputStream(trainingImagesFile.toString)
    )
  )

  val trainingLabelsFile = Paths.get("data/mnist/", "train-labels-idx1-ubyte.gz")
  val trainingLabels = new DataInputStream(
    new GZIPInputStream(
      new FileInputStream(trainingLabelsFile.toString)
    )
  )

  def readImage(
    rows:   Int,
    cols:   Int,
    stream: DataInputStream): DenseMatrix[Double] = {

    val m = new DenseMatrix[Double](rows, cols)

    for {
      i <- 0 until rows
      j <- 0 until cols
    } m(i, j) = stream.readUnsignedByte()

    m
  }

  def readImageFile(
    stream: DataInputStream,
    checknum: Int = 2051) = {

    val magicNumber = stream.readInt()
    if (magicNumber != checknum) {
      Left(throw new Exception("Wrong file mate"))
    } else {

      val total = stream.readInt()
      val rows = stream.readInt()
      val cols = stream.readInt()
      val images = Stream.continually(readImage(rows, cols, stream)).
        take(total)

      Right(images)
    }
  }

  def readLabelFile(
    stream: DataInputStream,
    checknum: Int = 2049) = {

    val magicNumber = stream.readInt()
    if (magicNumber != checknum) {
      Left(throw new Exception("Wrong file mate"))
    } else {

      val total = stream.readInt()
      val rows = stream.readInt()
      val cols = stream.readInt()
      val labels = Stream.continually(stream.readByte().toInt).
        take(total)

      Right(labels)
    }
  }

  val training = for {
    images <- readImageFile(trainingImages)
    labels <- readLabelFile(trainingLabels)
  } yield images zip labels
}

object ClassifyMnist extends App with ReadMnist {
//  val fitted = Classify.fit(ys, ks, tol, 10)
}




