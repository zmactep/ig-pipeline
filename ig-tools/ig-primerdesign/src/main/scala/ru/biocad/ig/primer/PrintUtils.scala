package ru.biocad.ig.primer
import Predef.{augmentString => _, _}

/**
 * Created by Kos on 22.02.14.
 */
object PrintUtils {
  def printSeqWithOffset(seq: List[Set[String]], offset: Int) {
    val prefix = List.fill(offset)(" ").mkString("")
    for (i <- 0 until 4) {
      val data = seq.map{s => if (i < s.size) s.toList(i) else " "}.mkString("")
      if (data.count(_ == ' ') < data.length) println(s"$prefix$data")
    }
  }

  def printOverlaps(strand: List[Set[String]], pos: List[Int], overlapSize: Int) = {
    for ((from, until) <- (0 :: pos).zip(pos)) {
      printSeqWithOffset(strand.slice(from, until + overlapSize), from)
    }
  }
}
