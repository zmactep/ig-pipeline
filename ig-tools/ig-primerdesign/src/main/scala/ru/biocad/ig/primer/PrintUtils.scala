package ru.biocad.ig.primer
import Predef.{augmentString => _, _}

/**
 * @author kfeodorov
 * @since 22.02.14
 */
object PrintUtils {
  def printSeqWithOffset(seq: List[Set[String]], offset: Int) {
    val prefix = List.fill(offset)(" ").mkString("")
    for (i <- 0 until 4) {
      val data = seq.map{s => if (i < s.size) s.toList(i) else " "}.mkString("")
      if (data.count(_ == ' ') < data.length) println(s"$prefix$data${if (i == 0) s" gc = ${DnaUtils.getGC(Option(data)).get}" else ""}")
    }
  }

  def printOverlaps(strand: Option[List[Set[String]]], pos: Option[List[Int]], overlapSize: Int) =
    strand.map{s => pos.map{ind: List[Int] => (0 :: ind).zip{ind :+ s.size}.foreach{case (from: Int, until: Int) => printSeqWithOffset(s.slice(from, until + overlapSize), from)}}}
}
