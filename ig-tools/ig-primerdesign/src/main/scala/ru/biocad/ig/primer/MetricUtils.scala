package ru.biocad.ig.primer

import scala.collection.mutable.ArrayBuffer

/**
 * @author kfeodorov
 * @since 21.03.14
 */
object MetricUtils {

  /**
   * 'Parrot' score - weight of longest common substring. C/G nucleotide has score of 3, A/T - of 2.
   * @param s1 sense strand
   * @param s2 reverse complement strand (or any other string)
   * @return weight of longest common substring + common substr
   */
  def similarityScore(s1: Option[String], s2: Option[String]): Option[(Int, String)] = {
    if (!s1.isDefined || !s2.isDefined) return None
    def matchScore(lhs: Char, rhs: Char): Int = if (lhs.toUpper != rhs.toUpper) 0
    else {if (List('C', 'G').contains(lhs.toUpper)) 3 else 2}

    val m = ArrayBuffer[Vector[Int]]()
    var prevRow: Vector[Int] = s1.get.map{c: Char => matchScore(c, s2.get.head)}.toVector
    m += prevRow

    for (c <- s2.get.drop(1)) {
      val currRow: Vector[Int] = s1.get.zipWithIndex.map{pair =>
        val score = matchScore(pair._1, c)
        if (score > 0) {
          score + m.last.lift(pair._2 - 1).getOrElse(0)
        } else score
      }.toVector
      prevRow = currRow
      m += currRow
    }

//    m.foreach(v => println(v.mkString("\t")))
//    println()

    val ans = m.flatten.max
    //find common substr given max score in matrix
    m.foreach{v: Vector[Int] =>
      val pos = v.indexOf(ans)
      if (pos >= 0) {
        var totalScore = ans
        val scoreMap = Map('C' -> 3, 'G' -> 3, 'A' -> 2, 'T' -> 2, 'U' -> 2)
        val substr = s1.get.reverse.drop(s1.get.size - pos - 1).takeWhile{c: Char => {if (totalScore > 0) {totalScore -= scoreMap(c); true} else false}}.reverse
        return Option((ans, substr))
      }
    }
    None
  }


  def hairpinScore(s: Option[String]): Option[Double] = similarityScore(s, DnaUtils.reverseComplementDNA(s)).flatMap{x => similarityScore(s, s) map {y => if (y == x) 0.5 else x._1.toDouble / y._1}}
  /**
   *
   * @param s
   * @param threshold
   * @return false if is None
   */
  def isHairpin(s: Option[String], threshold: Double): Boolean = hairpinScore(s) match {
    case Some(score) => score < threshold
    case None => false
  }

  /**
   * Calculates overlap score: similarityScore() is multiplied by one of coefs.
   * @param s1 1st primer
   * @param s2 2nd primer
   * @param coef - penalty. 1st coef is for most probable case where common substr is on the opposite ends:
   *             -----------
   *                     ------------
   *               2nd coef is for less common and undesired case:
   *              \___________/
   *               ___________
   *              /           \
   *              or:
   *              ___________/
   *              ___________
   *                         \
   * @return similarityScore() * coef
   */
  def overlapScore(s1: Option[String], s2: Option[String], coef: (Double, Double)): Option[Double] = similarityScore(s1, s2) match {
    case Some((score, substr)) =>
      val i1 = s1.get.indexOf(substr)
      val i2 = s2.get.indexOf(substr)
      val m1 = s1.get.size / 2
      val m2 = s2.get.size / 2
      //TODO take substr len into account
      Option(score * {
        if (i1 < m1 && i2 > m2 || i1 > m1 && i2 < m2) coef._1 else coef._2
      })
    case None => None
  }

  /**
   * Takes size/2 pieces from sense/antisense strands and compares it with the rest of strand.
   * Ensures that size/2 sticky end is unique across the strand
   * @param strand strand to parse
   * @param size primer size
   * @param coef penalty coefs for overlapScore()
   * @return
   */
  def primersOverlapScore(strand: Option[String], size: Int, coef: (Double, Double)): Option[Double] = strand map { s =>
    val rc = DnaUtils.reverseComplementRNA(strand).get

    def score(i: Int): Double = List(rc.take(i), rc.drop(i + size/2), s.take(i), s.drop(i + size/2)).filter(_.nonEmpty).
      map{ str => overlapScore(Option(str), Option(s.slice(i, i + size/2)), coef).get +
                  overlapScore(Option(str), Option(rc.slice(i, i + size/2)), coef).get}.sum
    (for {i <- 0 until s.length by size/2} yield score(i)).sum
  }
}
