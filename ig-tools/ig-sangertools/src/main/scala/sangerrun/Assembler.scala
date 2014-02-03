package sangerrun

import common.Algo

object Assembler {
  def getContig(reads: String*) : String = {
    reads.length match {
      case 0 => throw new IllegalArgumentException("Cannot assemble zero reads.")
      case 1 =>
        reads.head
      case 2 =>
        assemblePair(reads.head, reads.tail.head)
      case _ =>
        assembleAny(reads)
    }
  }

  private def assemblePair(seq1 : String, seq2 : String) : String = {
    val (pos1, pos2, len) = Algo.longestSubstr(seq1, seq2)
    seq1.substring(0, pos1 + len) + seq2.substring(pos2 + len)
  }

  private def assembleAny(reads: Seq[String]): String =
    throw new NotImplementedError()
}
