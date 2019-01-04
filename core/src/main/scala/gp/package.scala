package com.github.jonnylaw

import cats._
import breeze.stats.distributions.Rand

package object gp {
  implicit val randMonad = new Monad[Rand] {
    def pure[A](x: A): Rand[A] = Rand.always(x)
    def flatMap[A, B](fa: Rand[A])(f: A => Rand[B]): Rand[B] =
      fa flatMap f

    def tailRecM[A, B](a: A)(f: A => Rand[Either[A, B]]): Rand[B] =
      f(a).draw match {
        case Left(a1) => tailRecM(a1)(f)
        case Right(b) => Rand.always(b)
      }
  }
}
